"""Experiment 1: SPECIFICITY.
Score every molecule (blockers + non-blockers) by mean top-5 ECFP4 Tanimoto
to the hERG BLOCKER set, strict leave-one-out (exclude same InChIKey first).
ROC-AUC for separating blockers from non-blockers. Plus external drug sanity."""
import json, os, statistics
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs, RDLogger
RDLogger.DisableLog('rdApp.*')

HERE = os.path.dirname(os.path.abspath(__file__))
sets = json.load(open(os.path.join(HERE, "herg_sets.json")))
blockers = sets["blockers"]
nonblockers = sets["nonblockers"]

def featurize(dct):
    """returns list of (mid, inchikey, fp)"""
    out = []
    for mid, smi in dct.items():
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        try:
            ik = Chem.MolToInchiKey(m)
        except Exception:
            ik = mid
        fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048)
        out.append((mid, ik, fp))
    return out

print("Featurizing blockers...", file=__import__('sys').stderr)
B = featurize(blockers)
print("Featurizing non-blockers...", file=__import__('sys').stderr)
N = featurize(nonblockers)
print(f"valid blockers={len(B)} nonblockers={len(N)}", file=__import__('sys').stderr)

# reference blocker fps + parallel inchikey list
ref_fps = [fp for (_, _, fp) in B]
ref_iks = [ik for (_, ik, fp) in B]

def class_score(qfp, qik):
    """mean of top-5 Tanimoto to blocker set, excluding any ref with same InChIKey."""
    sims = DataStructs.BulkTanimotoSimilarity(qfp, ref_fps)
    # exclude self / duplicate structures
    filt = [s for s, ik in zip(sims, ref_iks) if ik != qik]
    filt.sort(reverse=True)
    top5 = filt[:5]
    return sum(top5) / len(top5)

blk_scores = []
for (mid, ik, fp) in B:
    blk_scores.append(class_score(fp, ik))
non_scores = []
for (mid, ik, fp) in N:
    non_scores.append(class_score(fp, ik))

def quartiles(v):
    v = sorted(v)
    n = len(v)
    def q(p):
        idx = p * (n - 1)
        lo = int(idx); hi = min(lo + 1, n - 1)
        return v[lo] + (v[hi] - v[lo]) * (idx - lo)
    return q(0.25), q(0.5), q(0.75)

def roc_auc(pos, neg):
    """rank-based AUC = P(score_pos > score_neg). Ties count 0.5."""
    allv = [(s, 1) for s in pos] + [(s, 0) for s in neg]
    allv.sort(key=lambda x: x[0])
    # assign ranks with tie handling (average ranks)
    ranks = [0.0] * len(allv)
    i = 0
    while i < len(allv):
        j = i
        while j < len(allv) and allv[j][0] == allv[i][0]:
            j += 1
        avg = (i + 1 + j) / 2.0  # average rank (1-based)
        for k in range(i, j):
            ranks[k] = avg
        i = j
    sum_pos = sum(r for r, (s, lab) in zip(ranks, allv) if lab == 1)
    npos = len(pos); nneg = len(neg)
    auc = (sum_pos - npos * (npos + 1) / 2.0) / (npos * nneg)
    return auc

auc = roc_auc(blk_scores, non_scores)
bq = quartiles(blk_scores)
nq = quartiles(non_scores)

# ---- External sanity: common non-cardiac drugs ----
EXTERNAL = {
 "metformin":"CN(C)C(=N)N=C(N)N",
 "aspirin":"CC(=O)OC1=CC=CC=C1C(=O)O",
 "acetaminophen":"CC(=O)NC1=CC=C(C=C1)O",
 "amoxicillin":"CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=C(C=C3)O)N)C(=O)O)C",
 "omeprazole":"CC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=CC=C(C=C3N2)OC",
 "ibuprofen":"CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
 "caffeine":"CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
 "lisinopril":"C1CC(N(C1)C(=O)C(CCCCN)NC(CCC2=CC=CC=C2)C(=O)O)C(=O)O",
 "metronidazole":"CC1=NC=C(N1CCO)[N+](=O)[O-]",
 "ranitidine":"CNC(=CNS(=O)(=O)C)NCCSCC1=CC=C(O1)CN(C)C",
 "amoxicillin2":"CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=C(C=C3)O)N)C(=O)O)C",
 "naproxen":"CC(C1=CC2=CC=C(C=C2C=C1)OC)C(=O)O",
 "gabapentin":"C1CCC(CC1)(CC(=O)O)CN",
 "hydrochlorothiazide":"C1=CC2=C(C=C1S(=O)(=O)N)S(=O)(=O)NCN2",
 "furosemide":"C1=CC=C(C(=C1)COC2=CC(=C(C=C2)S(=O)(=O)N)Cl)C(=O)O",
 "ciprofloxacin":"C1CC1N2C=C(C(=O)C3=CC(=C(C=C32)N4CCNCC4)F)C(=O)O",
 "levothyroxine":"C1=C(C=C(C(=C1I)OC2=CC(=C(C(=C2)I)O)I)I)CC(C(=O)O)N",
 "prednisone":"CC12CC(=O)C3C(C1CCC2(C(=O)CO)O)CCC4=CC(=O)C=CC34C",
 "warfarin":"CC(=O)CC(C1=CC=CC=C1)C2=C(C3=CC=CC=C3OC2=O)O",
 "ascorbic_acid":"C(C(C1C(=C(C(=O)O1)O)O)O)O",
}
ext_scores = {}
for n, smi in EXTERNAL.items():
    if n == "amoxicillin2":
        continue
    m = Chem.MolFromSmiles(smi)
    if m is None:
        ext_scores[n] = None; continue
    ik = Chem.MolToInchiKey(m)
    fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048)
    ext_scores[n] = round(class_score(fp, ik), 3)

# fraction of external that score below blocker median
bmed = bq[1]
ext_below = sum(1 for v in ext_scores.values() if v is not None and v < bmed)

result = {
    "counts": {"blockers_valid": len(B), "nonblockers_valid": len(N)},
    "auc": auc,
    "blocker_score_q25_med_q75": bq,
    "nonblocker_score_q25_med_q75": nq,
    "blocker_mean": statistics.mean(blk_scores),
    "nonblocker_mean": statistics.mean(non_scores),
    "external_drug_scores": ext_scores,
    "blocker_median": bmed,
    "external_below_blocker_median": f"{ext_below}/{len(ext_scores)}",
}
json.dump(result, open(os.path.join(HERE, "exp1_results.json"), "w"), indent=1)

print("\n===== EXPERIMENT 1: SPECIFICITY =====")
print(f"blockers scored: {len(B)}   non-blockers scored: {len(N)}")
print(f"ROC-AUC (blocker vs non-blocker) = {auc:.3f}")
print(f"Blocker class-score  q25/median/q75 = {bq[0]:.3f} / {bq[1]:.3f} / {bq[2]:.3f}")
print(f"Nonblk  class-score  q25/median/q75 = {nq[0]:.3f} / {nq[1]:.3f} / {nq[2]:.3f}")
print(f"\nExternal non-cardiac drugs (score to blocker class), median blocker={bmed:.3f}:")
for n, v in sorted(ext_scores.items(), key=lambda x: (x[1] is None, x[1])):
    flag = "" if (v is not None and v < bmed) else "  <-- HIGH"
    print(f"   {n:20s} {v}{flag}")
print(f"External below blocker median: {ext_below}/{len(ext_scores)}")
