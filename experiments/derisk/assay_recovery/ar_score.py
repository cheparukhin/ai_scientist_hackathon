import json, os, statistics
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

HERE = os.path.dirname(os.path.abspath(__file__))
data = json.load(open(os.path.join(HERE, "ar_actives.json")))

# ----- assay names (Layer-4 mapping) per target key -----
ASSAY = {
    "hERG_CHEMBL240":    "hERG patch-clamp",
    "Nav15_CHEMBL1980":  "Nav1.5 electrophysiology",
    "Cav12_CHEMBL1940":  "Cav1.2 electrophysiology",
    "5HT2B_CHEMBL1833":  "5-HT2B counter-screen",
    "5HT2A_CHEMBL224":   "5-HT2A binding",
    "D2_CHEMBL217":      "D2 binding",
    "M2_CHEMBL211":      "M2 muscarinic binding",
    "alpha1A_CHEMBL229": "alpha1A binding",
    "muOpioid_CHEMBL233":"mu-opioid binding",
    "SERT_CHEMBL228":    "SERT uptake",
    "AChE_CHEMBL220":    "AChE inhibition",
    "MAOA_CHEMBL1951":   "MAO-A inhibition",
    "COX1_CHEMBL221":    "COX-1 inhibition",
    "H1_CHEMBL231":      "H1 binding",
}

# ----- Baseline B: fixed, candidate-INDEPENDENT "standard practice" front-loading order -----
# Rationale: cardiac ion channels are universally front-loaded (hERG first, always);
# then core aminergic GPCR safety receptors (CNS/CV); then transporters/enzymes/specialized.
# 5-HT2B is placed among the aminergic GPCRs (Bowes-2012 core panel) - a CONSERVATIVE
# (baseline-favouring) placement, not buried at the bottom.
BASELINE_B_ORDER = [
    "hERG_CHEMBL240",     # 1 cardiac ion channel - always first
    "Nav15_CHEMBL1980",   # 2 cardiac
    "Cav12_CHEMBL1940",   # 3 cardiac
    "5HT2A_CHEMBL224",    # 4 aminergic GPCR
    "5HT2B_CHEMBL1833",   # 5 aminergic GPCR (valve counter-screen, Bowes core)
    "D2_CHEMBL217",       # 6
    "alpha1A_CHEMBL229",  # 7
    "M2_CHEMBL211",       # 8
    "H1_CHEMBL231",       # 9
    "muOpioid_CHEMBL233", # 10
    "SERT_CHEMBL228",     # 11
    "MAOA_CHEMBL1951",    # 12
    "AChE_CHEMBL220",     # 13
    "COX1_CHEMBL221",     # 14
]

QUERIES = {
    "terfenadine":  "CC(C)(C)C1=CC=C(C=C1)C(CCCN2CCC(CC2)C(C3=CC=CC=C3)(C4=CC=CC=C4)O)O",
    "cisapride":    "CO[C@H]1CN(CC[C@H]1NC(=O)C2=CC(=C(C=C2OC)N)Cl)CCCOC3=CC=C(C=C3)F",
    "astemizole":   "COC1=CC=C(C=C1)CCN2CCC(CC2)NC3=NC4=CC=CC=C4N3CC5=CC=C(C=C5)F",
    "thioridazine": "CN1CCCCC1CCN2C3=CC=CC=C3SC4=C2C=C(C=C4)SC",
    "pergolide":    "CCCN1C[C@@H](C[C@H]2[C@H]1CC3=CNC4=CC=CC2=C34)CSC",
    "fenfluramine": "CCNC(C)CC1=CC(=CC=C1)C(F)(F)F",
}
# killer target per query
TRUE = {
    "terfenadine": "hERG_CHEMBL240", "cisapride": "hERG_CHEMBL240",
    "astemizole": "hERG_CHEMBL240", "thioridazine": "hERG_CHEMBL240",
    "pergolide": "5HT2B_CHEMBL1833", "fenfluramine": "5HT2B_CHEMBL1833",
}
HERG_QUERIES = {"terfenadine", "cisapride", "astemizole", "thioridazine"}
HT2B_QUERIES = {"pergolide", "fenfluramine"}

# same 24-drug background as score.py
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
    "gabapentin":"C1CCC(CC1)(CC(=O)O)CN",
    "sertraline":"CNC1CCC(C2=CC=CC=C12)C3=CC(=C(C=C3)Cl)Cl",
    "ciprofloxacin":"C1CC1N2C=C(C(=O)C3=CC(=C(C=C32)N4CCNCC4)F)C(=O)O",
    "prednisone":"CC12CC(=O)C3C(C1CCC2(C(=O)CO)O)CCC4=CC(=O)C=CC34C",
    "levothyroxine":"C1=C(C=C(C(=C1I)OC2=CC(=C(C(=C2)I)O)I)I)CC(C(=O)O)N",
    "diazepam":"CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3",
    "ranitidine":"CNC(=CNS(=O)(=O)C)NCCSCC1=CC=C(O1)CN(C)C",
    "loratadine":"CCOC(=O)N1CCC(=C2C3=C(CCC4=C2N=CC=C4)C=C(C=C3)Cl)CC1",
    "naproxen":"CC(C1=CC2=CC=C(C=C2C=C1)OC)C(=O)O",
    "clopidogrel":"COC(=O)C(C1=CC=CC=C1Cl)N2CCC3=C(C2)C=CS3",
}

def fp(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None, None
    ik = Chem.MolToInchiKey(m)
    f = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048)
    return f, ik

print("Building target fingerprints...")
TARGET_FPS = {}
for tname, mols in data.items():
    lst = []
    for mid, smi in mols.items():
        f, ik = fp(smi)
        if f is not None:
            lst.append((ik, f))
    TARGET_FPS[tname] = lst
    print(f"  {tname}: {len(lst)} valid fps")

targets = list(data.keys())

def class_score(qfp, tname, exclude_iks):
    sims = []
    for ik, f in TARGET_FPS[tname]:
        if ik in exclude_iks:
            continue
        sims.append(DataStructs.TanimotoSimilarity(qfp, f))
    sims.sort(reverse=True)
    top5 = sims[:5]
    return statistics.mean(top5) if top5 else 0.0

QIK = {n: fp(s)[1] for n, s in QUERIES.items()}

def loo_exclude(qname):
    ex = set()
    ex.add(QIK[qname])
    grp = HERG_QUERIES if qname in HERG_QUERIES else HT2B_QUERIES
    for g in grp:
        ex.add(QIK[g])
    return ex

# background stats per target
BG = {t: [] for t in targets}
for n, s in BACKGROUND.items():
    f, ik = fp(s)
    if f is None:
        continue
    for t in targets:
        BG[t].append(class_score(f, t, set([ik])))
bg_stats = {t: (statistics.mean(BG[t]), statistics.pstdev(BG[t])) for t in targets}

def zscore(v, t):
    m, sd = bg_stats[t]
    return (v - m) / sd if sd > 0 else float('nan')

# per-query z per target
rows = {}
for qn, qs in QUERIES.items():
    qf, _ = fp(qs)
    ex = loo_exclude(qn)
    z = {}
    for t in targets:
        z[t] = zscore(class_score(qf, t, ex), t)
    rows[qn] = z

# ----- rankings -----
alpha_order = sorted(targets, key=lambda t: ASSAY[t].lower())   # Baseline A

def rank_of(ordered_list, killer):
    return ordered_list.index(killer) + 1

results = {"targets": targets, "assay": ASSAY,
           "baseline_A_alpha_order": [ASSAY[t] for t in alpha_order],
           "baseline_B_order": [ASSAY[t] for t in BASELINE_B_ORDER],
           "bg_stats": bg_stats, "z": rows, "per_drug": {}}

print("\n=== KILLER-ASSAY RANK per drug (of 14) ===")
print(f"{'drug':13} {'killer assay':26} {'OURS':>5} {'baseA':>6} {'baseB':>6}")
ours_top1 = ours_top3 = 0
bA_top1 = bA_top3 = 0
bB_top1 = bB_top3 = 0
for qn in QUERIES:
    killer = TRUE[qn]
    ours_order = sorted(targets, key=lambda t: rows[qn][t], reverse=True)
    r_ours = rank_of(ours_order, killer)
    r_A = rank_of(alpha_order, killer)
    r_B = rank_of(BASELINE_B_ORDER, killer)
    results["per_drug"][qn] = {
        "killer_target": killer, "killer_assay": ASSAY[killer],
        "rank_ours": r_ours, "rank_baselineA_alpha": r_A, "rank_baselineB_default": r_B,
        "killer_z": rows[qn][killer],
        "ours_ranked_assays": [ASSAY[t] for t in ours_order],
        "ours_top3_assays": [ASSAY[t] for t in ours_order[:3]],
    }
    ours_top1 += r_ours == 1; ours_top3 += r_ours <= 3
    bA_top1 += r_A == 1;      bA_top3 += r_A <= 3
    bB_top1 += r_B == 1;      bB_top3 += r_B <= 3
    print(f"{qn:13} {ASSAY[killer]:26} {r_ours:>5} {r_A:>6} {r_B:>6}")

n = len(QUERIES)
results["recovery"] = {
    "ours":      {"top1": ours_top1, "top3": ours_top3, "n": n},
    "baselineA": {"top1": bA_top1,  "top3": bA_top3,  "n": n},
    "baselineB": {"top1": bB_top1,  "top3": bB_top3,  "n": n},
}
print(f"\nTop-1 recovery: ours {ours_top1}/{n}, baseA {bA_top1}/{n}, baseB {bB_top1}/{n}")
print(f"Top-3 recovery: ours {ours_top3}/{n}, baseA {bA_top3}/{n}, baseB {bB_top3}/{n}")

# also report full z-vector per drug for transparency
print("\n=== killer-target z-score per drug ===")
for qn in QUERIES:
    print(f"  {qn:13} z(killer={ASSAY[TRUE[qn]]}) = {rows[qn][TRUE[qn]]:+.2f}")

json.dump(results, open(os.path.join(HERE, "results.json"), "w"), indent=1)
print("\nWrote results.json")
