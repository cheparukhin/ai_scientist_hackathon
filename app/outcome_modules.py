"""outcome_modules.py - M2 outcome breadth (LOWER evidence tier than the M1 off-target core).

Two model-predicted modules the validated M1 assays-to-culprit engine is blind to:
  - reactive_alerts(smiles): reactive-metabolite structural ALERTS (hypotheses, not
    predictions) - RDKit BRENK/NIH FilterCatalog + curated bioactivation SMARTS.
  - outcome_scores(smiles): hepatotox / mito read-across - SAME R4 backbone as core.py
    (mean-top5 ECFP4 Tanimoto to the endpoint's curated actives, z vs the SAME 24-drug
    background as core). This is similarity-ENRICHMENT, NOT P(harm).

HONESTY (BUILD2 §2): everything here is z-score / similarity, never a probability. Tag
tier="model-predicted"; render below and visually separated from the M1 plan; never merge
into the assays-to-culprit headline.
"""
import json
import os
import statistics

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit import RDLogger
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")

NBITS = 2048
FLAG_Z = 2.0          # a read-across endpoint is "flagged" when z >= this (mirrors core)
TIER = "model-predicted"

# BRENK descriptors that fire on STABLE amidine/guanidine/imine C=N (e.g. metformin's
# biguanide) - these are generic medicinal-chemistry "unwanted" filters, NOT
# reactive-metabolite bioactivation alerts, so we do not surface them here. Keeping them
# would over-call clean molecules; the curated bioactivation SMARTS remain the primary net.
SUPPRESS_CATALOG = {"imine_1", "imine_2"}


def _fp(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=NBITS)


class _OutcomeEngine:
    """Loads curated outcome actives + the shared 24-drug background once."""

    def __init__(self):
        self.actives = json.load(open(os.path.join(DATA, "outcome_actives.json")))
        self.alerts = json.load(open(os.path.join(DATA, "reactive_alerts.json")))
        background = json.load(open(os.path.join(DATA, "background.json")))

        # per-endpoint active fingerprints (+ keep meta for nearest-analog citations)
        self.endpoint_fps = {}
        for ep, drugs in self.actives.items():
            lst = []
            for d in drugs:
                f = _fp(d["smiles"])
                if f is not None:
                    lst.append((d["name"], f, d["citation"]))
            self.endpoint_fps[ep] = lst

        # SAME 24-drug background as core.py
        self.bg_fps = []
        for name, smi in background.items():
            if name.startswith("_"):
                continue
            f = _fp(smi)
            if f is not None:
                self.bg_fps.append(f)

        # background stats per endpoint (each bg drug's mean-top5 class score)
        self.bg_stats = {}
        for ep in self.actives:
            vals = [self._class_score(f, ep) for f in self.bg_fps]
            self.bg_stats[ep] = (statistics.mean(vals), statistics.pstdev(vals))

        # curated SMARTS (pre-compiled) + RDKit BRENK/NIH catalog
        self.curated = []
        for a in self.alerts:
            patt = Chem.MolFromSmarts(a["smarts"])
            if patt is not None:
                self.curated.append((a["name"], patt, a["note"], a["citation"]))
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.NIH)
        self.catalog = FilterCatalog(params)

    def _class_score(self, qfp, endpoint):
        fps = [f for _n, f, _c in self.endpoint_fps[endpoint]]
        if not fps:
            return 0.0
        sims = DataStructs.BulkTanimotoSimilarity(qfp, fps)
        sims.sort(reverse=True)
        top5 = sims[:5]
        return statistics.mean(top5) if top5 else 0.0

    def _z(self, raw, endpoint):
        m, sd = self.bg_stats[endpoint]
        return (raw - m) / sd if sd > 0 else 0.0

    def _nearest(self, qfp, endpoint):
        best = (None, 0.0, None)
        for name, f, cit in self.endpoint_fps[endpoint]:
            s = DataStructs.TanimotoSimilarity(qfp, f)
            if s > best[1]:
                best = (name, s, cit)
        return best

    def reactive(self, smiles):
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            return []
        hits = {}  # name(lower) -> dict, dedupe by name
        # curated bioactivation SMARTS
        for name, patt, note, cit in self.curated:
            if m.HasSubstructMatch(patt):
                hits.setdefault(name.lower(), {"name": name, "note": note, "citation": cit})
        # RDKit BRENK/NIH catalog
        for entry in self.catalog.GetMatches(m):
            nm = entry.GetDescription()
            if nm in SUPPRESS_CATALOG:
                continue
            key = nm.lower()
            if key not in hits:
                hits[key] = {
                    "name": nm,
                    "note": "structural alert (RDKit BRENK/NIH catalog) - potential "
                            "reactive-metabolite / bioactivation liability.",
                    "citation": "RDKit FilterCatalog (BRENK; NIH). Hypothesis - confirm experimentally.",
                }
        return list(hits.values())

    def scores(self, smiles):
        f = _fp(smiles)
        if f is None:
            return {}
        out = {}
        for ep in self.actives:
            raw = self._class_score(f, ep)
            z = self._z(raw, ep)
            nm, sim, cit = self._nearest(f, ep)
            out[ep] = {
                "z": round(z, 2),
                "raw": round(raw, 3),
                "nearest": {"name": nm, "sim": round(sim, 3), "citation": cit},
                "tier": TIER,
                "flagged": z >= FLAG_Z,
            }
        return out


_ENGINE = None


def get_outcome_engine():
    global _ENGINE
    if _ENGINE is None:
        _ENGINE = _OutcomeEngine()
    return _ENGINE


def reactive_alerts(smiles):
    """Reactive-metabolite structural ALERTS (not predictions). [{name, note, citation}]."""
    return get_outcome_engine().reactive(smiles)


def outcome_scores(smiles):
    """Read-across enrichment per endpoint. {endpoint: {z, raw, nearest, tier, flagged}}."""
    return get_outcome_engine().scores(smiles)


def outcome_panel(smiles):
    """Full M2 model-predicted panel (lower tier than M1). All z-based, no probabilities."""
    return {
        "reactive": reactive_alerts(smiles),
        "endpoints": outcome_scores(smiles),
        "tier": TIER,
        "caveat": "Model-predicted, lower evidence tier than the M1 off-target core "
                  "(less validated). Similarity-enrichment / structural alerts, not a "
                  "probability of harm; confirm experimentally.",
    }
