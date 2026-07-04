"""Experiment 3: MULTI-TARGET SPECIFICITY.
Replicate the hERG class-membership specificity test on 3 non-cardiac targets,
one per functional class, plus an optional 4th fallback.

Method (identical to hERG exp1):
  - Fetch ChEMBL activities with a pchembl_value (IC50/Ki/Kd/Potency).
  - Label by MAX pChEMBL per molecule: active >=6, inactive <=4.5, drop 4.5-6.
  - Dedupe by molecule; canonical SMILES; further dedupe structures by InChIKey.
  - Score every molecule by mean top-5 ECFP4 (Morgan r=2, 2048-bit) Tanimoto to
    the ACTIVE set, strict leave-one-out (remove same-InChIKey refs first).
  - ROC-AUC (actives vs inactives) via Mann-Whitney rank; score medians (IQR).
  - Sanity: score ~10 unrelated drugs, confirm they land in inactive range.

stdlib urllib only; AUC computed by hand.
"""
import json, urllib.request, urllib.parse, time, sys, os, statistics
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs, RDLogger
RDLogger.DisableLog('rdApp.*')

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = "https://www.ebi.ac.uk/chembl/api/data/activity.json"

TARGETS = [
    ("5-HT2B (aminergic GPCR)", "CHEMBL1833"),
    ("SERT/SLC6A4 (transporter)", "CHEMBL228"),
    ("AChE (enzyme)", "CHEMBL220"),
    ("MAO-A (enzyme, fallback)", "CHEMBL1951"),
]

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
 "ciprofloxacin":"C1CC1N2C=C(C(=O)C3=CC(=C(C=C32)N4CCNCC4)F)C(=O)O",
 "hydrochlorothiazide":"C1=CC2=C(C=C1S(=O)(=O)N)S(=O)(=O)NCN2",
 "gabapentin":"C1CCC(CC1)(CC(=O)O)CN",
}


def fetch_all(tid):
    per_mol = {}   # mid -> {"smiles":..., "pchembls":[...]}
    offset = 0
    limit = 1000
    pages = 0
    total = None
    while True:
        params = {
            "target_chembl_id": tid,
            "pchembl_value__isnull": "false",
            "limit": limit,
            "offset": offset,
        }
        url = BASE + "?" + urllib.parse.urlencode(params)
        for attempt in range(5):
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "feas-test"})
                with urllib.request.urlopen(req, timeout=90) as r:
                    data = json.load(r)
                break
            except Exception as e:
                print("  err", e, file=sys.stderr); time.sleep(3)
        else:
            break
        acts = data.get("activities", [])
        if not acts:
            break
        for a in acts:
            st = a.get("standard_type")
            if st not in ("IC50", "Ki", "Kd", "Potency"):
                continue
            smi = a.get("canonical_smiles")
            mid = a.get("molecule_chembl_id")
            pv = a.get("pchembl_value")
            if not (smi and mid and pv):
                continue
            try:
                pv = float(pv)
            except ValueError:
                continue
            d = per_mol.setdefault(mid, {"smiles": smi, "pchembls": []})
            d["pchembls"].append(pv)
        pm = data.get("page_meta", {})
        total = pm.get("total_count")
        offset += limit
        pages += 1
        print(f"  [{tid}] page {pages}: uniq mols {len(per_mol)}, total {total}", file=sys.stderr)
        if not pm.get("next"):
            break
        time.sleep(0.3)
    return per_mol


def label(per_mol):
    actives, inactives = {}, {}
    ambiguous = 0
    for mid, d in per_mol.items():
        mx = max(d["pchembls"])
        if mx >= 6.0:
            actives[mid] = d["smiles"]
        elif mx <= 4.5:
            inactives[mid] = d["smiles"]
        else:
            ambiguous += 1
    return actives, inactives, ambiguous


def featurize(dct):
    """returns list of (mid, inchikey, fp), deduped by InChIKey (keep first)."""
    out = []
    seen = set()
    for mid, smi in dct.items():
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        try:
            ik = Chem.MolToInchiKey(m)
        except Exception:
            ik = mid
        if ik in seen:
            continue
        seen.add(ik)
        fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048)
        out.append((mid, ik, fp))
    return out


def quartiles(v):
    v = sorted(v)
    n = len(v)
    def q(p):
        idx = p * (n - 1)
        lo = int(idx); hi = min(lo + 1, n - 1)
        return v[lo] + (v[hi] - v[lo]) * (idx - lo)
    return q(0.25), q(0.5), q(0.75)


def roc_auc(pos, neg):
    """rank-based AUC = P(score_pos > score_neg). Ties count 0.5 (avg ranks)."""
    allv = [(s, 1) for s in pos] + [(s, 0) for s in neg]
    allv.sort(key=lambda x: x[0])
    ranks = [0.0] * len(allv)
    i = 0
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
    if npos == 0 or nneg == 0:
        return None
    return (sum_pos - npos * (npos + 1) / 2.0) / (npos * nneg)


def run_target(name, tid):
    print(f"\n### Fetching {name} [{tid}] ...", file=sys.stderr)
    per_mol = fetch_all(tid)
    actives, inactives, ambiguous = label(per_mol)
    print(f"### Featurizing {name} ...", file=sys.stderr)
    A = featurize(actives)
    N = featurize(inactives)
    ref_fps = [fp for (_, _, fp) in A]
    ref_iks = [ik for (_, ik, fp) in A]

    def class_score(qfp, qik):
        sims = DataStructs.BulkTanimotoSimilarity(qfp, ref_fps)
        filt = [s for s, ik in zip(sims, ref_iks) if ik != qik]
        if not filt:
            return 0.0
        filt.sort(reverse=True)
        top5 = filt[:5]
        return sum(top5) / len(top5)

    act_scores = [class_score(fp, ik) for (mid, ik, fp) in A]
    inact_scores = [class_score(fp, ik) for (mid, ik, fp) in N]

    auc = roc_auc(act_scores, inact_scores) if (act_scores and inact_scores) else None
    aq = quartiles(act_scores) if act_scores else (None, None, None)
    iq = quartiles(inact_scores) if inact_scores else (None, None, None)

    # external drug sanity (scored to this target's active class)
    ext_scores = {}
    for n, smi in EXTERNAL.items():
        m = Chem.MolFromSmiles(smi)
        if m is None or not ref_fps:
            ext_scores[n] = None; continue
        ik = Chem.MolToInchiKey(m)
        efp = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048)
        ext_scores[n] = round(class_score(efp, ik), 3)
    amed = aq[1]
    ext_below = sum(1 for v in ext_scores.values() if v is not None and amed is not None and v < amed)

    res = {
        "name": name, "target_chembl_id": tid,
        "counts": {
            "unique_molecules": len(per_mol),
            "actives_labeled": len(actives),
            "inactives_labeled": len(inactives),
            "ambiguous_4.5_to_6": ambiguous,
            "actives_valid_dedup": len(A),
            "inactives_valid_dedup": len(N),
        },
        "auc": auc,
        "active_score_q25_med_q75": aq,
        "inactive_score_q25_med_q75": iq,
        "active_median": aq[1],
        "inactive_median": iq[1],
        "active_mean": statistics.mean(act_scores) if act_scores else None,
        "inactive_mean": statistics.mean(inact_scores) if inact_scores else None,
        "external_drug_scores": ext_scores,
        "external_below_active_median": f"{ext_below}/{len(ext_scores)}",
    }
    return res


def main():
    results = []
    # Always run the first 3; run MAO-A only if any of the 3 has too few inactives.
    core = TARGETS[:3]
    fallback = TARGETS[3]
    MIN_INACT = 15
    for name, tid in core:
        results.append(run_target(name, tid))
    need_fallback = any(r["counts"]["inactives_valid_dedup"] < MIN_INACT for r in results)
    if need_fallback:
        print("\n### A core target had too few inactives; adding MAO-A fallback.", file=sys.stderr)
        results.append(run_target(*fallback))

    out = {"targets": results, "herg_reference_auc": 0.894}
    json.dump(out, open(os.path.join(HERE, "exp3_results.json"), "w"), indent=1)

    print("\n\n===== EXPERIMENT 3: MULTI-TARGET SPECIFICITY =====")
    print(f"(hERG reference ROC-AUC = 0.894)\n")
    for r in results:
        c = r["counts"]
        aq = r["active_score_q25_med_q75"]; iq = r["inactive_score_q25_med_q75"]
        print(f"--- {r['name']} [{r['target_chembl_id']}] ---")
        print(f"  unique mols={c['unique_molecules']}  actives(>=6)={c['actives_valid_dedup']}  "
              f"inactives(<=4.5)={c['inactives_valid_dedup']}  ambiguous(dropped)={c['ambiguous_4.5_to_6']}")
        if r["auc"] is not None:
            print(f"  ROC-AUC (active vs inactive) = {r['auc']:.3f}")
            print(f"  active   score q25/med/q75 = {aq[0]:.3f} / {aq[1]:.3f} / {aq[2]:.3f}")
            print(f"  inactive score q25/med/q75 = {iq[0]:.3f} / {iq[1]:.3f} / {iq[2]:.3f}")
            print(f"  external drugs below active median: {r['external_below_active_median']}")
        else:
            print("  ROC-AUC = N/A (insufficient class data)")
        print()


if __name__ == "__main__":
    main()
