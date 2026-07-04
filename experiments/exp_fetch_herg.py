"""Fetch hERG (CHEMBL240) activities with pchembl_value, dedupe by molecule,
classify blockers (max pChEMBL>=6) vs non-blockers (max pChEMBL<=4.5).
Uses max pChEMBL per molecule so a compound potent in ANY assay is never
labelled a non-blocker."""
import json, urllib.request, urllib.parse, time, sys, os

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
TID = "CHEMBL240"

def fetch_all():
    per_mol = {}   # mid -> {"smiles":..., "pchembls":[...]}
    offset = 0
    limit = 1000
    pages = 0
    while True:
        params = {
            "target_chembl_id": TID,
            "pchembl_value__isnull": "false",
            "limit": limit,
            "offset": offset,
        }
        url = BASE + "?" + urllib.parse.urlencode(params)
        for attempt in range(4):
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
        offset += limit
        pages += 1
        print(f"  page {pages}: uniq mols {len(per_mol)}, total {pm.get('total_count')}", file=sys.stderr)
        if not pm.get("next"):
            break
        time.sleep(0.3)
    return per_mol

per_mol = fetch_all()

blockers = {}
nonblockers = {}
ambiguous = 0
for mid, d in per_mol.items():
    mx = max(d["pchembls"])
    if mx >= 6.0:
        blockers[mid] = d["smiles"]
    elif mx <= 4.5:
        nonblockers[mid] = d["smiles"]
    else:
        ambiguous += 1

out = {
    "blockers": blockers,
    "nonblockers": nonblockers,
    "counts": {
        "unique_molecules": len(per_mol),
        "blockers": len(blockers),
        "nonblockers": len(nonblockers),
        "ambiguous_4.5_to_6": ambiguous,
    },
}
json.dump(out, open(os.path.join(HERE, "herg_sets.json"), "w"))
print("DONE", out["counts"])
