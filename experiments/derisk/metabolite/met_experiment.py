"""
De-risk probe: does metabolite-aware scoring rescue the fenfluramine 5-HT2B miss?

Reuses score.py's exact scoring: mean top-5 Tanimoto (ECFP4, Morgan r=2, 2048 bit)
to a target's ChEMBL active class, reported as a z-score vs the same 24-drug
BACKGROUND used in score.py, with leave-one-out removal of the query InChIKey and
known mechanistic partners.

Outputs (this folder only):
  met_actives.json   ChEMBL cache (5-HT2B, hERG, SERT actives)
  results.json       all scores
"""
import json, os, statistics, urllib.request, urllib.parse, time, sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "met_actives.json")
BASE = "https://www.ebi.ac.uk/chembl/api/data/activity.json"

TARGETS = {
    "5HT2B_CHEMBL1833": "CHEMBL1833",
    "hERG_CHEMBL240":   "CHEMBL240",
    "SERT_CHEMBL228":   "CHEMBL228",
}

# ---- exact BACKGROUND from experiments/score.py (24 unrelated drugs) ----
BACKGROUND = {
 "metformin":"CN(C)C(=N)N=C(N)N","aspirin":"CC(=O)OC1=CC=CC=C1C(=O)O",
 "acetaminophen":"CC(=O)NC1=CC=C(C=C1)O",
 "atorvastatin":"CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4",
 "omeprazole":"CC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=CC=C(C=C3N2)OC",
 "amoxicillin":"CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=C(C=C3)O)N)C(=O)O)C",
 "ibuprofen":"CC(C)CC1=CC=C(C=C1)C(C)C(=O)O","caffeine":"CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
 "lisinopril":"C1CC(N(C1)C(=O)C(CCCCN)NC(CCC2=CC=CC=C2)C(=O)O)C(=O)O",
 "warfarin":"CC(=O)CC(C1=CC=CC=C1)C2=C(C3=CC=CC=C3OC2=O)O",
 "furosemide":"C1=CC=C(C(=C1)COC2=CC(=C(C=C2)S(=O)(=O)N)Cl)C(=O)O",
 "hydrochlorothiazide":"C1=CC2=C(C=C1S(=O)(=O)N)S(=O)(=O)NCN2",
 "metoprolol":"CC(C)NCC(COC1=CC=C(C=C1)CCOC)O",
 "amlodipine":"CCOC(=O)C1=C(NC(=C(C1C2=CC=CC=C2Cl)C(=O)OC)C)COCCN",
 "simvastatin":"CCC(C)(C)C(=O)OC1CC(C=C2C1C(C(C=C2)C)CCC3CC(CC(=O)O3)O)C",
 "gabapentin":"C1CCC(CC1)(CC(=O)O)CN","sertraline":"CNC1CCC(C2=CC=CC=C12)C3=CC(=C(C=C3)Cl)Cl",
 "ciprofloxacin":"C1CC1N2C=C(C(=O)C3=CC(=C(C=C32)N4CCNCC4)F)C(=O)O",
 "prednisone":"CC12CC(=O)C3C(C1CCC2(C(=O)CO)O)CCC4=CC(=O)C=CC34C",
 "levothyroxine":"C1=C(C=C(C(=C1I)OC2=CC(=C(C(=C2)I)O)I)I)CC(C(=O)O)N",
 "diazepam":"CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3",
 "ranitidine":"CNC(=CNS(=O)(=O)C)NCCSCC1=CC=C(O1)CN(C)C",
 "loratadine":"CCOC(=O)N1CCC(=C2C3=C(CCC4=C2N=CC=C4)C=C(C=C3)Cl)CC1",
 "naproxen":"CC(C1=CC2=CC=C(C=C2C=C1)OC)C(=O)O",
 "clopidogrel":"COC(=O)C(C1=CC=CC=C1Cl)N2CCC3=C(C2)C=CS3",
}

# ---- parent / metabolite test pairs ----
# label -> (smiles, target, direction-note)
MOLS = {
 # fenfluramine: parent inactive at 5-HT2B, metabolite norfenfluramine is the agonist
 "fenfluramine":    ("CCNC(C)CC1=CC(=CC=C1)C(F)(F)F",           "5HT2B_CHEMBL1833"),
 "norfenfluramine": ("CC(Cc1cccc(c1)C(F)(F)F)N",                "5HT2B_CHEMBL1833"),
 # benfluorex (Mediator): prodrug -> norfenfluramine; parent inactive, metab is the 5-HT2B agonist
 "benfluorex":      ("O=C(OCCNC(C)Cc1cccc(c1)C(F)(F)F)c1ccccc1","5HT2B_CHEMBL1833"),
 # counter-example: terfenadine is the hERG blocker; its active metabolite fexofenadine is NON-cardiotoxic
 "terfenadine":     ("CC(C)(C)C1=CC=C(C=C1)C(CCCN2CCC(CC2)C(C3=CC=CC=C3)(C4=CC=CC=C4)O)O", "hERG_CHEMBL240"),
 "fexofenadine":    ("OC(c1ccccc1)(c1ccccc1)C1CCN(CCCC(O)c2ccc(cc2)C(C)(C)C(=O)O)CC1", "hERG_CHEMBL240"),
}

# metabolite relationships (parent -> metabolite) for ensemble scoring
PAIRS = {
 "fenfluramine": "norfenfluramine",
 "benfluorex":   "norfenfluramine",
 "terfenadine":  "fexofenadine",
}

def fetch_target(tid, want=900):
    mols = {}
    offset = 0; limit = 1000; pages = 0
    while len(mols) < want and pages < 8:
        params = {"target_chembl_id": tid, "pchembl_value__gte": 6,
                  "limit": limit, "offset": offset}
        url = BASE + "?" + urllib.parse.urlencode(params)
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "derisk-metabolite"})
            with urllib.request.urlopen(req, timeout=60) as r:
                data = json.load(r)
        except Exception as e:
            print("  err", e, file=sys.stderr); time.sleep(3); pages += 1; continue
        acts = data.get("activities", [])
        if not acts: break
        for a in acts:
            if a.get("standard_type") not in ("IC50", "Ki", "Kd"): continue
            smi = a.get("canonical_smiles"); mid = a.get("molecule_chembl_id")
            if smi and mid and mid not in mols: mols[mid] = smi
        pm = data.get("page_meta", {})
        offset += limit; pages += 1
        print(f"  {tid}: page {pages}, uniq {len(mols)}, total_acts {pm.get('total_count')}", file=sys.stderr)
        if not pm.get("next"): break
        time.sleep(0.5)
    return mols

# ---- load or fetch ChEMBL actives ----
if os.path.exists(CACHE):
    data = json.load(open(CACHE))
    print("Loaded cache", {k: len(v) for k, v in data.items()}, file=sys.stderr)
else:
    data = {}
    for name, tid in TARGETS.items():
        print("Fetching", name, tid, file=sys.stderr)
        data[name] = fetch_target(tid)
        print(name, "->", len(data[name]), file=sys.stderr)
    json.dump(data, open(CACHE, "w"))

def fp(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None: return None, None
    ik = Chem.MolToInchiKey(m)
    f = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048)
    return f, ik

# ---- build target fingerprints ----
TARGET_FPS = {}
for tname, mols in data.items():
    lst = []
    for mid, smi in mols.items():
        f, ik = fp(smi)
        if f is not None: lst.append((ik, f))
    TARGET_FPS[tname] = lst
    print(f"  {tname}: {len(lst)} valid fps", file=sys.stderr)

targets = list(data.keys())

def class_score(qfp, tname, exclude_iks):
    sims = []
    for ik, f in TARGET_FPS[tname]:
        if ik in exclude_iks: continue
        sims.append(DataStructs.TanimotoSimilarity(qfp, f))
    sims.sort(reverse=True)
    return statistics.mean(sims[:5]), (max(sims) if sims else 0.0), len(sims)

# ---- background stats per target (exclude the bg drug's own InChIKey, as score.py does) ----
BG = {t: [] for t in targets}
for n, s in BACKGROUND.items():
    f, ik = fp(s)
    if f is None: continue
    for t in targets:
        m5, mx, nc = class_score(f, t, {ik})
        BG[t].append((n, m5, mx))
bg_stats = {t: (statistics.mean([x[1] for x in BG[t]]),
                statistics.pstdev([x[1] for x in BG[t]]),
                max([x[1] for x in BG[t]])) for t in targets}

def z(v, t):
    m, sd, _ = bg_stats[t]
    return (v - m) / sd if sd > 0 else float('nan')

# ---- query InChIKeys, for LOO exclusion ----
QIK = {name: fp(smi)[1] for name, (smi, _t) in MOLS.items()}

# mechanistic partners to co-remove during LOO, per target class (mirrors score.py intent)
# 5-HT2B: all our 5-HT2B queries + pergolide; hERG: all our hERG queries + score.py's hERG set
HT2B_PARTNERS = {"fenfluramine", "norfenfluramine", "benfluorex"}
HERG_PARTNERS = {"terfenadine", "fexofenadine"}
# named external partners not in MOLS but present in ChEMBL classes -> exclude by InChIKey
EXTRA_PARTNER_SMILES = {
 "5HT2B_CHEMBL1833": {"pergolide": "CCCN1C[C@@H](C[C@H]2[C@H]1CC3=CNC4=CC=CC2=C34)CSC",
                      "cabergoline": "CCNC(=O)N(CCCN(C)C)C(=O)[C@H]1CN([C@@H]2Cc3c[nH]c4cccc(c34)C2=C1)CC=C"},
 "hERG_CHEMBL240": {"cisapride": "CO[C@H]1CN(CC[C@H]1NC(=O)C2=CC(=C(C=C2OC)N)Cl)CCCOC3=CC=C(C=C3)F",
                    "astemizole": "COC1=CC=C(C=C1)CCN2CCC(CC2)NC3=NC4=CC=CC=C4N3CC5=CC=C(C=C5)F"},
}
EXTRA_IK = {t: {nm: fp(s)[1] for nm, s in d.items()} for t, d in EXTRA_PARTNER_SMILES.items()}

def loo_exclude(label, target):
    ex = set()
    ex.add(QIK[label])
    partners = HT2B_PARTNERS if target == "5HT2B_CHEMBL1833" else (HERG_PARTNERS if target == "hERG_CHEMBL240" else set())
    for p in partners:
        if p in QIK: ex.add(QIK[p])
    for nm, ik in EXTRA_IK.get(target, {}).items():
        if ik: ex.add(ik)
    return ex

# ---- score every molecule on its target (and also 5-HT2B/hERG for context) ----
scored = {}
for label, (smi, target) in MOLS.items():
    qf, qik = fp(smi)
    row = {}
    for t in targets:
        ex = loo_exclude(label, t)
        # is the query's own InChIKey in this target class? (self-match check)
        present = any(ik == qik for ik, _ in TARGET_FPS[t])
        m5, mx, nc = class_score(qf, t, ex)
        row[t] = {"mean_top5": m5, "max": mx, "n_compared": nc,
                  "z": z(m5, t), "self_in_class": present, "n_excluded": len(ex)}
    scored[label] = {"smiles": smi, "inchikey": qik, "primary_target": target, "scores": row}

# ---- ensemble: parent vs metabolite on the parent's PRIMARY target ----
ensembles = {}
for parent, metab in PAIRS.items():
    tgt = MOLS[parent][1]
    pz = scored[parent]["scores"][tgt]["z"]
    mz = scored[metab]["scores"][tgt]["z"]
    ensembles[parent] = {
        "target": tgt, "metabolite": metab,
        "parent_z": pz, "metabolite_z": mz,
        "ensemble_max_z": max(pz, mz),
        "ensemble_mean_z": (pz + mz) / 2,
        "parent_mean_top5": scored[parent]["scores"][tgt]["mean_top5"],
        "metabolite_mean_top5": scored[metab]["scores"][tgt]["mean_top5"],
    }

out = {
    "targets": targets,
    "n_actives": {t: len(TARGET_FPS[t]) for t in targets},
    "bg_stats": {t: {"mean": bg_stats[t][0], "sd": bg_stats[t][1], "max": bg_stats[t][2]} for t in targets},
    "scored": scored,
    "ensembles": ensembles,
}
json.dump(out, open(os.path.join(HERE, "results.json"), "w"), indent=1)

# ---- report ----
print("\n=== BACKGROUND (24 drugs) mean-top5 class score ===")
for t in targets:
    m, sd, mx = bg_stats[t]
    print(f"  {t}: mean={m:.3f} sd={sd:.3f} max={mx:.3f}  (n_actives={len(TARGET_FPS[t])})")

print("\n=== PER-MOLECULE z on PRIMARY target (LOO applied) ===")
for label, (smi, tgt) in MOLS.items():
    sc = scored[label]["scores"][tgt]
    print(f"  {label:16s} -> {tgt:16s}  mean_top5={sc['mean_top5']:.3f}  z={sc['z']:+.2f}  "
          f"self_in_class={sc['self_in_class']}  n_excl={sc['n_excluded']}")

print("\n=== PARENT vs METABOLITE vs ENSEMBLE (z on parent's primary target) ===")
print(f"{'parent':14s}{'target':18s}{'parent_z':>10s}{'metab_z':>10s}{'ens_max':>10s}{'ens_mean':>10s}")
for parent, e in ensembles.items():
    print(f"{parent:14s}{e['target']:18s}{e['parent_z']:>10.2f}{e['metabolite_z']:>10.2f}"
          f"{e['ensemble_max_z']:>10.2f}{e['ensemble_mean_z']:>10.2f}")

print("\nWrote results.json + met_actives.json")
