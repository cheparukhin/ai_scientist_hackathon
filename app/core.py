"""core.py - THE scoring engine (R4 backbone, reused from experiments/).

score_candidate(smiles, metabolite_smiles=None) -> dict

Pipeline (do NOT reinvent the math - this is the validated R4 backbone from
experiments/score.py / experiments/derisk/expanded/exp_score.py):
  1. Parse + canonicalize. Run the descriptor-box ABSTAIN gate
     (experiments/derisk/abstain/calibrate_ad.py). No nearest-neighbour rule.
  2. Per panel target: mean-top5 ECFP4 (Morgan r=2, 2048-bit) Tanimoto to that
     target's ChEMBL actives, z-scored vs the 24-drug background. **No leave-one-out**
     (a real novel candidate removes nothing).
  3. known-analog flag: max single Tanimoto to any reference_failures drug >= 0.5.
  4. Metabolite: if given, score parent AND metabolite, aggregate per target by MAX.

Honesty guardrails (BUILD.md sec 2) live in M2/UI copy; this module returns z-scores,
N-fold enrichment and ranks - never a probability of harm.
"""
import json
import math
import os
import statistics

from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, Descriptors
from rdkit import DataStructs
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")

# ---- abstain-gate thresholds (calibrated in calibrate_ad.py; 0/31 false-abstain,
#      0/11 false-accept on the test set). Do NOT loosen without re-validating. ----
MIN_HEAVY = 6
MW_MAX = 800.0
FRACHALO_MAX = 0.30
FRACO_MAX = 0.33
NBITS = 2048
Z_FLOOR = 2.0            # informational weak-coverage flag (NOT an abstain trigger)
KNOWN_ANALOG_T = 0.50    # known-analog flag threshold (max single Tanimoto to a ref failure)
FLAG_Z = 2.0             # a target is "flagged" (worth surfacing) when z >= this

METALS = set([3, 4, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 38,
              39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 55, 56, 57, 72, 73, 74, 75, 76,
              77, 78, 79, 80, 81, 82, 83])

# canonical example molecules for the demo (SMILES from experiments/derisk/drugs.json
# + calibrate_ad.py). "novel_kinase_like" is a drug-like molecule that is NOT a known
# failed drug - it exercises the novel-chemotype path the engine is built for.
EXAMPLES = {
    "rimonabant":       "Cc1c(C(=O)NN2CCCCC2)nn(-c2ccc(Cl)cc2Cl)c1-c1ccc(Cl)cc1",
    "pergolide":        "CCCN1C[C@H](CSC)C[C@@H]2c3cccc4[nH]cc(c34)C[C@H]21",
    "terfenadine":      "CC(C)(C)c1ccc(C(O)CCCN2CCC(C(O)(c3ccccc3)c3ccccc3)CC2)cc1",
    "cyclosporine":     "CCC1NC(=O)C(C(O)C(C)CC=CC)N(C)C(=O)C(C(C)C)N(C)C(=O)C(CC(C)C)N(C)C(=O)C(CC(C)C)N(C)C(=O)C(C)NC(=O)C(C)NC(=O)C(CC(C)C)N(C)C(=O)C(C(C)C)NC(=O)C(CC(C)C)N(C)C(=O)CN(C)C1=O",
    "novel_arylpiperazine": "O=C(Nc1ccc(N2CCN(c3ncccn3)CC2)cc1)c1ccc(Cl)cc1",
}


def _fp_and_desc(smi):
    """Return (mol, fp, descriptors dict) or (None, None, None) if unparseable."""
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None, None, None
    fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=NBITS)
    heavy = m.GetNumHeavyAtoms()
    nums = [a.GetAtomicNum() for a in m.GetAtoms() if a.GetAtomicNum() > 1]
    halo = sum(1 for z in nums if z in (9, 17, 35, 53))
    nO = sum(1 for z in nums if z == 8)
    desc = {
        "canonical_smiles": Chem.MolToSmiles(m),
        "heavy": heavy,
        "MW": round(Descriptors.MolWt(m), 1),
        "logP": round(Crippen.MolLogP(m), 2),
        "fracHalo": round(halo / heavy, 3) if heavy else 0.0,
        "fracO": round(nO / heavy, 3) if heavy else 0.0,
        "metal": any(a.GetAtomicNum() in METALS for a in m.GetAtoms()),
    }
    return m, fp, desc


def _abstain_check(desc):
    """Descriptor-box scope gate. Returns (abstain: bool, gate, reason)."""
    if desc["metal"]:
        return True, "A", "metal-containing - not a scoreable organic small molecule (out of scope)"
    if desc["heavy"] < MIN_HEAVY:
        return True, "A", f"<{MIN_HEAVY} heavy atoms - too small to score (out of scope)"
    if desc["MW"] > MW_MAX:
        return True, "B", f"MW {desc['MW']:.0f} > {MW_MAX:.0f} - physicochemically outside the reference set (large peptide / macrocycle)"
    if desc["fracHalo"] > FRACHALO_MAX:
        return True, "B", f"halogen fraction {desc['fracHalo']:.2f} > {FRACHALO_MAX} - outside the reference chemical space"
    if desc["fracO"] > FRACO_MAX:
        return True, "B", f"oxygen fraction {desc['fracO']:.2f} > {FRACO_MAX} - outside the reference chemical space (sugar / polyol / polyether)"
    return False, None, "in-domain"


class _Engine:
    """Loads panel actives / background / reference fingerprints once (heavy)."""

    def __init__(self):
        panel = json.load(open(os.path.join(DATA, "panel.json")))
        actives = json.load(open(os.path.join(DATA, "panel_actives.json")))
        background = json.load(open(os.path.join(DATA, "background.json")))
        refs = json.load(open(os.path.join(DATA, "reference_failures.json")))

        self.panel = panel
        self.targets = list(panel.keys())

        # per-target active fingerprints
        self.target_fps = {}
        for t in self.targets:
            fps = []
            for _mid, smi in actives.get(t, {}).items():
                m = Chem.MolFromSmiles(smi)
                if m is not None:
                    fps.append(AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=NBITS))
            self.target_fps[t] = fps

        # background stats per target (each bg drug's mean-top5 class score)
        bg = {t: [] for t in self.targets}
        for name, smi in background.items():
            if name.startswith("_"):
                continue
            m = Chem.MolFromSmiles(smi)
            if m is None:
                continue
            f = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=NBITS)
            for t in self.targets:
                bg[t].append(self._class_score(f, t))
        self.bg_stats = {t: (statistics.mean(v), statistics.pstdev(v)) for t, v in bg.items()}

        # reference-failure fingerprints (for known-analog flag + evidence)
        self.ref = []
        for name, r in refs.items():
            m = Chem.MolFromSmiles(r["smiles"])
            if m is None:
                continue
            f = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=NBITS)
            self.ref.append((name, f, r))

    def _class_score(self, qfp, target):
        fps = self.target_fps[target]
        if not fps:
            return 0.0
        sims = DataStructs.BulkTanimotoSimilarity(qfp, fps)
        sims.sort(reverse=True)
        top5 = sims[:5]
        return statistics.mean(top5) if top5 else 0.0

    def _z(self, raw, target):
        m, sd = self.bg_stats[target]
        return (raw - m) / sd if sd > 0 else 0.0

    def score_fp(self, fp):
        """Per-target raw mean-top5 and z for one fingerprint."""
        out = {}
        for t in self.targets:
            raw = self._class_score(fp, t)
            out[t] = {"raw": raw, "z": self._z(raw, t)}
        return out

    def nearest_ref(self, fp):
        """(name, sim, meta) of the most similar reference-failure drug."""
        best = (None, 0.0, None)
        for name, rfp, meta in self.ref:
            s = DataStructs.TanimotoSimilarity(fp, rfp)
            if s > best[1]:
                best = (name, s, meta)
        return best

    def refs_for_target(self, target_key):
        """All reference-failure drugs whose culprit target is this panel target."""
        return [meta for _n, _f, meta in self.ref if meta.get("culprit_target_key") == target_key]

    def max_ref_sim_at_target(self, fp, target_key):
        """Max single Tanimoto from fp to any reference failure whose culprit is target_key."""
        best = 0.0
        for _n, rfp, meta in self.ref:
            if meta.get("culprit_target_key") == target_key:
                best = max(best, DataStructs.TanimotoSimilarity(fp, rfp))
        return best


# ---- Layer-3b severity re-weight (consequence weighting; gentle so it never breaks
#      the validated z-driven assays-to-culprit ranking) ----
SEV_WEIGHT = {"high": 1.0, "med": 0.9, "low": 0.8}


def _priority(z, severity):
    """0-100 priority INDEX (NOT a probability). Saturating map of severity-weighted z."""
    eff = max(z, 0.0) * SEV_WEIGHT.get(severity, 0.9)
    return eff, int(round(100 * (1 - math.exp(-eff / 3.0))))


def build_plan(result, fp=None):
    """Turn a score_candidate() 'ok' result into a reordered assay plan with evidence.

    Returns None if the result abstained. Otherwise:
      {status, rows:[...ordered by our_rank...], headline:{...}, flagged:[keys]}
    Each row: target_key, assay_name, organ, severity, z, priority(int 0-100),
      our_rank, default_rank, delta, action(no-go|counter-screen|monitor),
      flagged(bool), grounded(bool), evidence:[failed-drug dicts], marginal_value, note.
    """
    if result.get("status") != "ok":
        return None

    eng = get_engine()
    panel = eng.panel
    baseline = json.load(open(os.path.join(DATA, "baseline_order.json")))["order"]
    default_rank = {t: i + 1 for i, t in enumerate(baseline)}

    if fp is None:
        _m, fp, _d = _fp_and_desc(result["canonical_smiles"])

    scored = result["targets"]

    # effective (severity-weighted) priority per target -> our ranking
    eff = {}
    for t in eng.targets:
        e, _p = _priority(scored[t]["z"], panel[t]["severity"])
        eff[t] = e
    our_order = sorted(eng.targets, key=lambda t: (eff[t], scored[t]["z"]), reverse=True)
    our_rank = {t: i + 1 for i, t in enumerate(our_order)}

    rows = []
    for t in our_order:
        z = scored[t]["z"]
        sev = panel[t]["severity"]
        _e, prio = _priority(z, sev)
        flagged = z >= FLAG_Z
        # grounded = candidate structurally resembles a KNOWN failed drug at THIS target
        grounded = (fp is not None) and (eng.max_ref_sim_at_target(fp, t) >= KNOWN_ANALOG_T)
        evidence = [
            {"name": m["name"], "organ": m["organ"], "tier": m["tier"], "citation": m["citation"]}
            for m in eng.refs_for_target(t)
        ]

        # action tag from severity + flagging + known-failure grounding
        if not flagged:
            action = "monitor"
        elif sev == "high":
            action = "no-go" if grounded else "counter-screen"
        elif sev == "med":
            action = "counter-screen"
        else:
            action = "monitor"

        rows.append({
            "target_key": t,
            "assay_name": panel[t]["assay_name"],
            "organ": panel[t]["organ"],
            "severity": sev,
            "marginal_value": panel[t]["marginal_value"],
            "note": panel[t]["note"],
            "z": round(z, 2),
            "priority": prio,
            "our_rank": our_rank[t],
            "default_rank": default_rank[t],
            "delta": default_rank[t] - our_rank[t],
            "action": action,
            "flagged": flagged,
            "grounded": grounded,
            "evidence": evidence,
        })

    flagged_keys = [r["target_key"] for r in rows if r["flagged"]]

    # headline = the top flagged assay (the one our reordering moves to catch the
    # program-ending liability first). Falls back to rank-1 assay if nothing flagged.
    head_row = next((r for r in rows if r["flagged"]), rows[0])
    headline = {
        "assay_name": head_row["assay_name"],
        "target_key": head_row["target_key"],
        "our_rank": head_row["our_rank"],
        "default_rank": head_row["default_rank"],
        "delta": head_row["delta"],
        "organ": head_row["organ"],
        "action": head_row["action"],
        "marginal_value": head_row["marginal_value"],
        "any_flagged": bool(flagged_keys),
    }

    return {"status": "ok", "rows": rows, "headline": headline, "flagged": flagged_keys}


_ENGINE = None


def get_engine():
    global _ENGINE
    if _ENGINE is None:
        _ENGINE = _Engine()
    return _ENGINE


def score_candidate(smiles, metabolite_smiles=None):
    """Score a candidate SMILES against the 18-target safety panel.

    Returns a dict:
      abstain:  {status:"abstain", reason, gate, descriptors}
      ok:       {status:"ok", canonical_smiles, descriptors, targets:{key:{raw,z}},
                 best_target, best_z, flags:{weak_coverage, known_analog},
                 nearest_analog:{name,sim,...}, metabolite:{...}|None}
    """
    m, fp, desc = _fp_and_desc(smiles)
    if m is None:
        return {"status": "abstain", "gate": "A", "reason": "unparseable SMILES", "input": smiles}

    abstain, gate, reason = _abstain_check(desc)
    if abstain:
        return {"status": "abstain", "gate": gate, "reason": reason,
                "descriptors": desc, "canonical_smiles": desc["canonical_smiles"],
                "input": smiles}

    eng = get_engine()
    parent = eng.score_fp(fp)

    # metabolite MAX-aggregation (rule from experiments/derisk/metabolite/)
    metabolite_info = None
    scored = parent
    if metabolite_smiles:
        mm, mfp, mdesc = _fp_and_desc(metabolite_smiles)
        if mm is not None:
            met = eng.score_fp(mfp)
            merged = {}
            for t in eng.targets:
                if met[t]["z"] > parent[t]["z"]:
                    merged[t] = dict(met[t]); merged[t]["from"] = "metabolite"
                else:
                    merged[t] = dict(parent[t]); merged[t]["from"] = "parent"
            scored = merged
            metabolite_info = {"canonical_smiles": mdesc["canonical_smiles"], "aggregation": "MAX(parent, metabolite) per target"}
        else:
            metabolite_info = {"error": "metabolite SMILES unparseable - scored parent only"}

    # known-analog flag: max single Tanimoto to any reference failure (parent, and
    # metabolite if provided) >= threshold
    na_name, na_sim, na_meta = eng.nearest_ref(fp)
    if metabolite_smiles:
        mm, mfp, _ = _fp_and_desc(metabolite_smiles)
        if mm is not None:
            n2, s2, m2 = eng.nearest_ref(mfp)
            if s2 > na_sim:
                na_name, na_sim, na_meta = n2, s2, m2
    known_analog = na_sim >= KNOWN_ANALOG_T

    best_target = max(eng.targets, key=lambda t: scored[t]["z"])
    best_z = scored[best_target]["z"]

    return {
        "status": "ok",
        "canonical_smiles": desc["canonical_smiles"],
        "descriptors": desc,
        "targets": scored,
        "best_target": best_target,
        "best_z": best_z,
        "flags": {
            "weak_coverage": best_z < Z_FLOOR,
            "known_analog": known_analog,
        },
        "nearest_analog": {
            "name": na_name,
            "sim": round(na_sim, 3),
            "culprit_target": na_meta.get("culprit_target") if na_meta else None,
            "organ": na_meta.get("organ") if na_meta else None,
            "citation": na_meta.get("citation") if na_meta else None,
        },
        "metabolite": metabolite_info,
    }
