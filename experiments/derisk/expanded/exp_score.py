import json, os, statistics
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

HERE = os.path.dirname(os.path.abspath(__file__))
data = json.load(open(os.path.join(HERE, "exp_actives.json")))
DRUGS = json.load(open(os.path.join(HERE, "..", "drugs.json")))

# ----- assay names (Layer-4 mapping) per target key -----
ASSAY = {
    "hERG_CHEMBL240":     "hERG patch-clamp",
    "Nav15_CHEMBL1980":   "Nav1.5 electrophysiology",
    "Cav12_CHEMBL1940":   "Cav1.2 electrophysiology",
    "5HT2A_CHEMBL224":    "5-HT2A binding",
    "5HT2B_CHEMBL1833":   "5-HT2B counter-screen",
    "5HT3_CHEMBL1899":    "5-HT3 binding",
    "D2_CHEMBL217":       "D2 binding",
    "M2_CHEMBL211":       "M2 muscarinic binding",
    "alpha1A_CHEMBL229":  "alpha1A binding",
    "H1_CHEMBL231":       "H1 binding",
    "muOpioid_CHEMBL233": "mu-opioid binding",
    "SERT_CHEMBL228":     "SERT uptake",
    "NET_CHEMBL222":      "NET uptake",
    "DAT_CHEMBL238":      "DAT uptake",
    "CB1_CHEMBL218":      "CB1 counter-screen",
    "AChE_CHEMBL220":     "AChE inhibition",
    "MAOA_CHEMBL1951":    "MAO-A inhibition",
    "COX1_CHEMBL221":     "COX-1 inhibition",
}

# map culprit_target string in drugs.json -> panel target key
CULPRIT_TO_KEY = {
    "5-HT2B": "5HT2B_CHEMBL1833",
    "CB1":    "CB1_CHEMBL218",
    "5-HT3":  "5HT3_CHEMBL1899",
    "NET":    "NET_CHEMBL222",
    "Cav1.2": "Cav12_CHEMBL1940",
    "hERG":   "hERG_CHEMBL240",
}

# ----- Baseline B: fixed, candidate-INDEPENDENT "standard practice" order (18) -----
# hERG, Nav1.5, Cav1.2 front-loaded (cardiac channels, as real panels do); buried
# off-targets sit low - honest.
BASELINE_B_ORDER = [
    "hERG_CHEMBL240",     # 1
    "Nav15_CHEMBL1980",   # 2
    "Cav12_CHEMBL1940",   # 3
    "5HT2A_CHEMBL224",    # 4
    "D2_CHEMBL217",       # 5
    "M2_CHEMBL211",       # 6
    "alpha1A_CHEMBL229",  # 7
    "H1_CHEMBL231",       # 8
    "muOpioid_CHEMBL233", # 9
    "SERT_CHEMBL228",     # 10
    "5HT2B_CHEMBL1833",   # 11
    "5HT3_CHEMBL1899",    # 12
    "NET_CHEMBL222",      # 13
    "DAT_CHEMBL238",      # 14
    "CB1_CHEMBL218",      # 15
    "AChE_CHEMBL220",     # 16
    "MAOA_CHEMBL1951",    # 17
    "COX1_CHEMBL221",     # 18
]

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

targets = list(ASSAY.keys())
n_actives = {t: len(TARGET_FPS[t]) for t in targets}

def class_score(qfp, tname, exclude_iks):
    sims = []
    for ik, f in TARGET_FPS[tname]:
        if ik in exclude_iks:
            continue
        sims.append(DataStructs.TanimotoSimilarity(qfp, f))
    sims.sort(reverse=True)
    top5 = sims[:5]
    return statistics.mean(top5) if top5 else 0.0

# query inchikeys (recompute from SMILES so exclusion matches ChEMBL InChIKeys)
QIK = {}
QFP = {}
for dname, d in DRUGS.items():
    f, ik = fp(d["smiles"])
    QFP[dname] = f
    QIK[dname] = ik

# group drugs by culprit_target for strict LOO (query + mechanistic partners)
by_culprit = {}
for dname, d in DRUGS.items():
    by_culprit.setdefault(d["culprit_target"], []).append(dname)

def loo_exclude(qname):
    ct = DRUGS[qname]["culprit_target"]
    ex = set()
    ex.add(QIK[qname])
    for partner in by_culprit[ct]:
        ex.add(QIK[partner])
    return ex

# ----- background stats per target (each bg drug excludes only itself) -----
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

# ----- per-query z + raw top5 mean per target -----
rows = {}      # z per target
raw_rows = {}  # raw top5-mean per target
for dname, d in DRUGS.items():
    qf = QFP[dname]
    ex = loo_exclude(dname)
    z = {}
    raw = {}
    for t in targets:
        cs = class_score(qf, t, ex)
        raw[t] = cs
        z[t] = zscore(cs, t)
    rows[dname] = z
    raw_rows[dname] = raw

def rank_of(ordered_list, killer):
    return ordered_list.index(killer) + 1

# ----- Baseline C: base-rate leave-one-out ranking (per query) -----
# rank assays by how often each target is the culprit across the OTHER 19 drugs.
# tie-break: Baseline B order (stable, defensible).
b_order_index = {t: i for i, t in enumerate(BASELINE_B_ORDER)}
def baseline_C_order(qname):
    counts = {t: 0 for t in targets}
    for other, d in DRUGS.items():
        if other == qname:
            continue
        key = CULPRIT_TO_KEY[d["culprit_target"]]
        counts[key] += 1
    return sorted(targets, key=lambda t: (-counts[t], b_order_index[t])), counts

results = {
    "targets": targets,
    "assay": ASSAY,
    "n_actives": n_actives,
    "baseline_B_order": [ASSAY[t] for t in BASELINE_B_ORDER],
    "bg_stats": bg_stats,
    "z": rows,
    "raw_top5mean": raw_rows,
    "per_drug": {},
}

cats = ["buried", "cardiac_herg"]
# recovery accumulators per category
rec = {c: {m: {"top1": 0, "top3": 0, "n": 0} for m in ["ours", "B", "C"]} for c in cats}
rank_sum = {c: {m: 0 for m in ["ours", "B", "C"]} for c in cats}
money = {c: 0 for c in cats}          # ours top3 AND B not top3
money_drugs = {c: [] for c in cats}

print("\n=== KILLER-ASSAY RANK per drug (of 18) ===")
print(f"{'drug':16} {'cat':13} {'killer assay':26} {'OURS':>5}{'B':>4}{'C':>4} {'z':>7}")
for dname, d in DRUGS.items():
    cat = d["category"]
    killer = CULPRIT_TO_KEY[d["culprit_target"]]
    ours_order = sorted(targets, key=lambda t: rows[dname][t], reverse=True)
    r_ours = rank_of(ours_order, killer)
    r_B = rank_of(BASELINE_B_ORDER, killer)
    c_order, c_counts = baseline_C_order(dname)
    r_C = rank_of(c_order, killer)
    kz = rows[dname][killer]

    results["per_drug"][dname] = {
        "category": cat,
        "culprit_target": d["culprit_target"],
        "killer_target": killer,
        "killer_assay": ASSAY[killer],
        "rank_ours": r_ours,
        "rank_B_default": r_B,
        "rank_C_baserate": r_C,
        "killer_z": kz,
        "killer_raw_top5mean": raw_rows[dname][killer],
        "ours_ranked_assays": [ASSAY[t] for t in ours_order],
        "ours_top3_assays": [ASSAY[t] for t in ours_order[:3]],
        "z_vector": {ASSAY[t]: rows[dname][t] for t in targets},
    }

    rec[cat]["ours"]["n"] += 1
    rec[cat]["B"]["n"] += 1
    rec[cat]["C"]["n"] += 1
    rec[cat]["ours"]["top1"] += r_ours == 1
    rec[cat]["ours"]["top3"] += r_ours <= 3
    rec[cat]["B"]["top1"] += r_B == 1
    rec[cat]["B"]["top3"] += r_B <= 3
    rec[cat]["C"]["top1"] += r_C == 1
    rec[cat]["C"]["top3"] += r_C <= 3
    rank_sum[cat]["ours"] += r_ours
    rank_sum[cat]["B"] += r_B
    rank_sum[cat]["C"] += r_C
    if r_ours <= 3 and r_B > 3:
        money[cat] += 1
        money_drugs[cat].append(dname)

    print(f"{dname:16} {cat:13} {ASSAY[killer]:26} {r_ours:>5}{r_B:>4}{r_C:>4} {kz:>+7.2f}")

# expected assays-to-culprit (mean killer rank) per category
mean_rank = {}
for c in cats:
    ncat = rec[c]["ours"]["n"]
    mean_rank[c] = {m: rank_sum[c][m] / ncat for m in ["ours", "B", "C"]}

results["recovery_by_category"] = rec
results["mean_expected_assays_to_culprit"] = mean_rank
results["money_number"] = money
results["money_drugs"] = money_drugs

print("\n=== RECOVERY by category (killer in top-K) ===")
for c in cats:
    ncat = rec[c]["ours"]["n"]
    print(f"\n[{c}] n={ncat}")
    print(f"  Top-1:  Ours {rec[c]['ours']['top1']}/{ncat}  B {rec[c]['B']['top1']}/{ncat}  C {rec[c]['C']['top1']}/{ncat}")
    print(f"  Top-3:  Ours {rec[c]['ours']['top3']}/{ncat}  B {rec[c]['B']['top3']}/{ncat}  C {rec[c]['C']['top3']}/{ncat}")
    print(f"  Mean expected assays-to-culprit:  Ours {mean_rank[c]['ours']:.2f}  B {mean_rank[c]['B']:.2f}  C {mean_rank[c]['C']:.2f}")

print("\n=== THE MONEY NUMBER (Ours top-3 AND Baseline B NOT top-3) ===")
for c in cats:
    print(f"  [{c}] {money[c]}/{rec[c]['ours']['n']}  -> {money_drugs[c]}")

json.dump(results, open(os.path.join(HERE, "results.json"), "w"), indent=1)
print("\nWrote results.json")
