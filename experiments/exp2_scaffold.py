"""Experiment 2: SCAFFOLD GENERALIZATION (leave-one-scaffold-cluster-out).
Cluster hERG blockers by Bemis-Murcko scaffold. For the top ~5 clusters,
hold out the ENTIRE cluster; score each held-out blocker using ONLY blockers
of a DIFFERENT scaffold. Compare held-out-blocker scores vs non-blocker
distribution, and compute AUC per holdout against the same reduced reference."""
import json, os, sys, statistics
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import DataStructs, RDLogger
RDLogger.DisableLog('rdApp.*')

HERE = os.path.dirname(os.path.abspath(__file__))
sets = json.load(open(os.path.join(HERE, "herg_sets.json")))
blockers = sets["blockers"]
nonblockers = sets["nonblockers"]

def featurize(dct, with_scaffold=False):
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
        scaf = None
        if with_scaffold:
            try:
                sm = MurckoScaffold.GetScaffoldForMol(m)
                scaf = Chem.MolToSmiles(sm) if sm is not None else ""
            except Exception:
                scaf = ""
        out.append((mid, ik, fp, scaf))
    return out

print("Featurizing...", file=sys.stderr)
B = featurize(blockers, with_scaffold=True)
N = featurize(nonblockers)
non_fps = [fp for (_, _, fp, _) in N]
non_iks = [ik for (_, ik, fp, _) in N]

# ---- cluster blockers by scaffold ----
scaf_groups = defaultdict(list)
for idx, (mid, ik, fp, scaf) in enumerate(B):
    key = scaf if scaf else "<none>"
    scaf_groups[key].append(idx)

# rank clusters by size (ignore empty-scaffold bucket which is acyclic mixed)
clusters = sorted(
    [(k, v) for k, v in scaf_groups.items() if k != "<none>" and k != ""],
    key=lambda kv: len(kv[1]), reverse=True)
top_clusters = clusters[:5]

print(f"Total distinct Murcko scaffolds among blockers: {len(scaf_groups)}", file=sys.stderr)

def top5_mean(qfp, qik, ref_fps, ref_iks):
    sims = DataStructs.BulkTanimotoSimilarity(qfp, ref_fps)
    filt = [s for s, ik in zip(sims, ref_iks) if ik != qik]
    filt.sort(reverse=True)
    top5 = filt[:5]
    return sum(top5) / len(top5)

def quartiles(v):
    v = sorted(v); n = len(v)
    def q(p):
        idx = p * (n - 1); lo = int(idx); hi = min(lo + 1, n - 1)
        return v[lo] + (v[hi] - v[lo]) * (idx - lo)
    return q(0.25), q(0.5), q(0.75)

def roc_auc(pos, neg):
    allv = [(s, 1) for s in pos] + [(s, 0) for s in neg]
    allv.sort(key=lambda x: x[0])
    ranks = [0.0] * len(allv); i = 0
    while i < len(allv):
        j = i
        while j < len(allv) and allv[j][0] == allv[i][0]:
            j += 1
        avg = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[k] = avg
        i = j
    sum_pos = sum(r for r, (s, lab) in zip(ranks, allv) if lab == 1)
    npos = len(pos); nneg = len(neg)
    return (sum_pos - npos * (npos + 1) / 2.0) / (npos * nneg)

# in-set (Exp1-style) scores for the held-out members, for degradation comparison
all_blk_fps = [fp for (_, _, fp, _) in B]
all_blk_iks = [ik for (_, ik, fp, _) in B]

# non-blocker scores vs FULL blocker set (baseline reference)
non_scores_full = [top5_mean(fp, ik, all_blk_fps, all_blk_iks) for (_, ik, fp, _) in N]
nq_full = quartiles(non_scores_full)

results = {"nonblocker_full_ref_q25_med_q75": nq_full, "clusters": []}
print("\n===== EXPERIMENT 2: SCAFFOLD LEAVE-ONE-CLUSTER-OUT =====")
print(f"Non-blocker score (vs full blocker ref) median = {nq_full[1]:.3f}")
print(f"(Exp1 blocker in-set median was ~0.664, AUC 0.894)\n")

for scaf, idxs in top_clusters:
    member_set = set(idxs)
    # reduced reference = blockers NOT in this cluster
    ref_fps = [B[i][2] for i in range(len(B)) if i not in member_set]
    ref_iks = [B[i][1] for i in range(len(B)) if i not in member_set]

    # held-out blockers scored vs reduced ref (novel-scaffold condition)
    held_scores = [top5_mean(B[i][2], B[i][1], ref_fps, ref_iks) for i in idxs]
    # same members' in-set score (full ref, LOO) for degradation delta
    inset_scores = [top5_mean(B[i][2], B[i][1], all_blk_fps, all_blk_iks) for i in idxs]
    # non-blockers scored vs SAME reduced ref for a fair AUC
    non_scores_red = [top5_mean(fp, ik, ref_fps, ref_iks) for (_, ik, fp, _) in N]

    auc_holdout = roc_auc(held_scores, non_scores_red)
    hq = quartiles(held_scores)
    iq = quartiles(inset_scores)
    frac_above_nonmed = sum(1 for s in held_scores if s > nq_full[1]) / len(held_scores)

    rec = {
        "scaffold_smiles": scaf,
        "cluster_size": len(idxs),
        "inset_score_median": iq[1],
        "holdout_score_q25_med_q75": hq,
        "holdout_auc_vs_nonblockers": auc_holdout,
        "median_drop_inset_to_holdout": iq[1] - hq[1],
        "frac_holdout_above_nonblocker_median": frac_above_nonmed,
    }
    results["clusters"].append(rec)
    print(f"scaffold: {scaf}")
    print(f"  size={len(idxs)}  in-set median={iq[1]:.3f}  ->  holdout median={hq[1]:.3f}  (drop {iq[1]-hq[1]:+.3f})")
    print(f"  holdout score q25/med/q75 = {hq[0]:.3f}/{hq[1]:.3f}/{hq[2]:.3f}")
    print(f"  holdout AUC vs non-blockers = {auc_holdout:.3f}")
    print(f"  frac held-out above non-blocker median = {frac_above_nonmed:.2f}\n")

# pooled: all held-out blockers across the 5 clusters vs their reduced-ref nonblocker scores
# (use full-ref nonblockers as common negative baseline for a single pooled AUC)
pooled_held = []
for scaf, idxs in top_clusters:
    member_set = set(idxs)
    ref_fps = [B[i][2] for i in range(len(B)) if i not in member_set]
    ref_iks = [B[i][1] for i in range(len(B)) if i not in member_set]
    pooled_held += [top5_mean(B[i][2], B[i][1], ref_fps, ref_iks) for i in idxs]
pooled_auc = roc_auc(pooled_held, non_scores_full)
results["pooled_holdout_auc_vs_full_nonblockers"] = pooled_auc
results["pooled_holdout_median"] = quartiles(pooled_held)[1]
print(f"POOLED (top-5 clusters held out) blocker median={quartiles(pooled_held)[1]:.3f}")
print(f"POOLED holdout AUC vs non-blockers = {pooled_auc:.3f}")

json.dump(results, open(os.path.join(HERE, "exp2_results.json"), "w"), indent=1)
