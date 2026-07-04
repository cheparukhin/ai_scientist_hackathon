"""ABLATION: is R4 (ChEMBL class-aggregation) NECESSARY, or does dumb 2D
nearest-failed-drug match it on the assay-ranking task?

Three methods rank the 18 panel assays for each of the 20 query drugs; we
record the KILLER-ASSAY RANK (rank of the assay mapped to the drug's
culprit_target) under each.

  R4          : mean-top5 ECFP4 Tanimoto to the whole ChEMBL active class,
                z-scored vs the 24-drug background (score.py engine). Strict LOO:
                query + all same-culprit_target drugs.json partners removed
                (skeleton-InChIKey match) from EVERY ChEMBL class before scoring.
  NAIVE-2D    : per target, MAX ECFP4 Tanimoto between the query and any OTHER
                drugs.json drug whose culprit_target is that target. Retrieval
                against the ~19 reference drugs only. Self excluded.
  NAIVE-2D+CH : per target, Tanimoto to ONE fixed random ChEMBL active of that
                target (isolates whether AGGREGATION over the class is what matters).

Tie handling for ranks: midrank (average rank) -- neutral, standard. NAIVE has
many exact-0 scores (targets with no reference drug), so its killer often lands
in a tie block; midrank neither favours nor punishes.
"""
import json, os, statistics, random
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

HERE = os.path.dirname(os.path.abspath(__file__))
data = json.load(open(os.path.join(HERE, "abl_actives.json")))
drugs = json.load(open(os.path.join(HERE, "..", "drugs.json")))

# ---- target key -> named assay (Layer-4 mapping); assay name == target name ----
ASSAY = {
    "hERG_CHEMBL240": "hERG", "Nav15_CHEMBL1980": "Nav1.5", "Cav12_CHEMBL1940": "Cav1.2",
    "5HT2A_CHEMBL224": "5-HT2A", "5HT2B_CHEMBL1833": "5-HT2B", "5HT3_CHEMBL1899": "5-HT3",
    "D2_CHEMBL217": "D2", "M2_CHEMBL211": "M2", "alpha1A_CHEMBL229": "alpha1A",
    "H1_CHEMBL231": "H1", "muOpioid_CHEMBL233": "mu-opioid", "SERT_CHEMBL228": "SERT",
    "NET_CHEMBL222": "NET", "DAT_CHEMBL238": "DAT", "CB1_CHEMBL218": "CB1",
    "AChE_CHEMBL220": "AChE", "MAOA_CHEMBL1951": "MAO-A", "COX1_CHEMBL221": "COX-1",
}
targets = list(data.keys())
assert set(targets) == set(ASSAY), (set(targets) ^ set(ASSAY))

# chembl id -> target key (to map each drug's culprit to a panel target)
CHEMBLID_TO_TKEY = {t.split("_")[1]: t for t in targets}

def fp(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None, None
    ik = Chem.MolToInchiKey(m)
    return AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048), ik

def skel(ik):
    return ik.split("-")[0] if ik else None   # 14-char connectivity block

# ---- build ChEMBL target fingerprints (store skeleton for LOO) ----
print("Building target fingerprints...")
TARGET_FPS = {}
for tname, mols in data.items():
    lst = []
    for mid, smi in mols.items():
        f, ik = fp(smi)
        if f is not None:
            lst.append((skel(ik), f))
    TARGET_FPS[tname] = lst
    print(f"  {tname}: {len(lst)} valid fps")

# ---- query drugs from drugs.json ----
Q = {}          # name -> {fp, ik, skel, culprit_tkey, category}
for name, d in drugs.items():
    f, ik = fp(d["smiles"])
    tkey = CHEMBLID_TO_TKEY.get(d["culprit_chembl_id"])
    assert tkey is not None, f"{name} culprit {d['culprit_chembl_id']} not in panel"
    Q[name] = {"fp": f, "ik": ik, "skel": skel(ik), "tkey": tkey,
               "category": d["category"]}

# same-culprit partner groups (drugs.json only)
BY_TKEY = {}
for name, q in Q.items():
    BY_TKEY.setdefault(q["tkey"], []).append(name)

# ---- 24-drug background (identical to score.py) ----
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

# ================= R4 engine =================
def class_score(qfp, tname, exclude_skels):
    sims = [DataStructs.TanimotoSimilarity(qfp, f)
            for sk, f in TARGET_FPS[tname] if sk not in exclude_skels]
    sims.sort(reverse=True)
    top5 = sims[:5]
    return statistics.mean(top5) if top5 else 0.0

# background stats per target (exclude only the bg drug's own skeleton)
BG = {t: [] for t in targets}
for n, s in BACKGROUND.items():
    f, ik = fp(s)
    if f is None:
        continue
    for t in targets:
        BG[t].append(class_score(f, t, {skel(ik)}))
bg_stats = {t: (statistics.mean(BG[t]), statistics.pstdev(BG[t])) for t in targets}

def zscore(v, t):
    m, sd = bg_stats[t]
    return (v - m) / sd if sd > 0 else float('nan')

def r4_loo_exclude(qname):
    """query + all same-culprit_target drugs.json partners, by skeleton."""
    q = Q[qname]
    ex = {q["skel"]}
    for partner in BY_TKEY[q["tkey"]]:
        ex.add(Q[partner]["skel"])
    return ex

# ================= NAIVE-2D =================
def naive2d_score(qname, tname, exclude_partners=False):
    """MAX ECFP4 Tanimoto to OTHER drugs.json drugs whose culprit == tname.
    exclude_partners=True => STRICT variant that also drops same-culprit
    partners of the query (degenerate: empties the killer target)."""
    q = Q[qname]
    refs = BY_TKEY.get(tname, [])
    sims = []
    for r in refs:
        if r == qname:
            continue
        if exclude_partners and Q[r]["tkey"] == q["tkey"]:
            continue
        sims.append(DataStructs.TanimotoSimilarity(q["fp"], Q[r]["fp"]))
    return max(sims) if sims else 0.0

# ================= NAIVE-2D + ChEMBL (single random active) =================
rng = random.Random(42)
ONE_ACTIVE = {t: rng.choice(TARGET_FPS[t]) for t in targets}  # (skel, fp)
def naive_chembl_score(qname, tname):
    q = Q[qname]
    sk, f = ONE_ACTIVE[tname]
    if sk == q["skel"]:        # never self-match
        # pick another
        for sk2, f2 in TARGET_FPS[tname]:
            if sk2 != q["skel"]:
                f = f2; break
    return DataStructs.TanimotoSimilarity(q["fp"], f)

# ---- ranking helper: midrank of the killer (descending score) ----
def killer_midrank(scores, killer):
    kv = scores[killer]
    greater = sum(1 for t, v in scores.items() if v > kv)
    equal = sum(1 for t, v in scores.items() if v == kv)   # includes killer
    return greater + (equal + 1) / 2.0                      # midrank

# ================= run all drugs =================
per_drug = {}
for qname in Q:
    q = Q[qname]
    killer = q["tkey"]
    ex = r4_loo_exclude(qname)
    r4 = {t: zscore(class_score(q["fp"], t, ex), t) for t in targets}
    nv = {t: naive2d_score(qname, t) for t in targets}                      # partners retained
    nv_strict = {t: naive2d_score(qname, t, exclude_partners=True) for t in targets}
    nc = {t: naive_chembl_score(qname, t) for t in targets}
    per_drug[qname] = {
        "category": q["category"],
        "killer_target": killer, "killer_assay": ASSAY[killer],
        "n_partners_retained": len([r for r in BY_TKEY[killer] if r != qname]),
        "rank_R4": killer_midrank(r4, killer),
        "rank_NAIVE": killer_midrank(nv, killer),
        "rank_NAIVE_strictLOO": killer_midrank(nv_strict, killer),
        "rank_NAIVE_chembl": killer_midrank(nc, killer),
        "R4_killer_z": r4[killer],
        "NAIVE_killer_maxsim": nv[killer],
        "NAIVE_chembl_killer_sim": nc[killer],
        "R4_top3": [ASSAY[t] for t in sorted(targets, key=lambda t: r4[t], reverse=True)[:3]],
        "NAIVE_top3": [ASSAY[t] for t in sorted(targets, key=lambda t: nv[t], reverse=True)[:3]],
    }

# ================= aggregate metrics =================
def agg(names, key):
    ranks = [per_drug[n][key] for n in names]
    top1 = sum(1 for r in ranks if r <= 1.0)     # midrank<=1 means strict rank 1, no tie above
    top3 = sum(1 for r in ranks if r <= 3.0)
    return {"top1": top1, "top3": top3, "mean_rank": round(statistics.mean(ranks), 2), "n": len(names)}

buried = [n for n in Q if Q[n]["category"] == "buried"]
cardiac = [n for n in Q if Q[n]["category"] == "cardiac_herg"]
alln = list(Q)

summary = {}
for split, names in [("all", alln), ("buried", buried), ("cardiac_herg", cardiac)]:
    summary[split] = {
        "R4": agg(names, "rank_R4"),
        "NAIVE_2D": agg(names, "rank_NAIVE"),
        "NAIVE_2D_strictLOO": agg(names, "rank_NAIVE_strictLOO"),
        "NAIVE_2D_chembl": agg(names, "rank_NAIVE_chembl"),
    }

# head-to-head verdict counts (R4 top3 & NAIVE not; and vice-versa) -- partner-retained NAIVE
r4_wins = [n for n in alln if per_drug[n]["rank_R4"] <= 3 and per_drug[n]["rank_NAIVE"] > 3]
naive_wins = [n for n in alln if per_drug[n]["rank_NAIVE"] <= 3 and per_drug[n]["rank_R4"] > 3]
naive_ties_or_beats = [n for n in alln if per_drug[n]["rank_NAIVE"] <= per_drug[n]["rank_R4"]]

verdict = {
    "R4_top3_NAIVE_not": r4_wins,
    "NAIVE_top3_R4_not": naive_wins,
    "NAIVE_ties_or_beats_R4": naive_ties_or_beats,
    "n_drugs": len(alln),
}

out = {"targets": targets, "assay": ASSAY, "bg_stats": bg_stats,
       "per_drug": per_drug, "summary": summary, "verdict": verdict,
       "n_actives": {t: len(TARGET_FPS[t]) for t in targets}}
json.dump(out, open(os.path.join(HERE, "results.json"), "w"), indent=1)

# ================= print report =================
print("\n=== PER-DRUG killer-assay rank (of 18) ===")
hdr = f"{'drug':16}{'cat':13}{'killer':9}{'#prt':>5}{'R4':>6}{'NAIVE':>7}{'delta':>7}{'R4z':>7}{'Nsim':>6}"
print(hdr)
for n in cardiac + buried:
    p = per_drug[n]
    d = p["rank_NAIVE"] - p["rank_R4"]
    print(f"{n:16}{p['category']:13}{p['killer_assay']:9}{p['n_partners_retained']:>5}"
          f"{p['rank_R4']:>6.1f}{p['rank_NAIVE']:>7.1f}{d:>+7.1f}{p['R4_killer_z']:>7.2f}{p['NAIVE_killer_maxsim']:>6.2f}")

print("\n=== RECOVERY (top1 / top3 / mean-rank) ===")
for split in ("all", "buried", "cardiac_herg"):
    s = summary[split]
    print(f"[{split}] n={s['R4']['n']}")
    for m in ("R4", "NAIVE_2D", "NAIVE_2D_strictLOO", "NAIVE_2D_chembl"):
        x = s[m]
        print(f"    {m:20} top1={x['top1']:>2}/{x['n']}  top3={x['top3']:>2}/{x['n']}  mean_rank={x['mean_rank']}")

print("\n=== VERDICT (partner-retained NAIVE) ===")
print(f"R4 top-3 while NAIVE-2D NOT: {len(r4_wins)}  {r4_wins}")
print(f"NAIVE-2D top-3 while R4 NOT: {len(naive_wins)}  {naive_wins}")
print(f"NAIVE ties-or-beats R4:      {len(naive_ties_or_beats)}  {naive_ties_or_beats}")
print("\nWrote results.json")
