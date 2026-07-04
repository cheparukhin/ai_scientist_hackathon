import json, urllib.request, urllib.parse, time, sys, os

OUT = os.path.dirname(os.path.abspath(__file__))
BASE = "https://www.ebi.ac.uk/chembl/api/data/activity.json"

TARGETS = {
    "hERG_CHEMBL240": "CHEMBL240",
    "5HT2B_CHEMBL1833": "CHEMBL1833",
    "BSEP_CHEMBL6020": "CHEMBL6020",
}

def fetch_target(tid, want=1000):
    # collect activities with pchembl>=6, standard_type in IC50/Ki/Kd
    mols = {}  # molecule_chembl_id -> canonical_smiles
    offset = 0
    limit = 1000
    pages = 0
    while len(mols) < want and pages < 8:
        params = {
            "target_chembl_id": tid,
            "pchembl_value__gte": 6,
            "limit": limit,
            "offset": offset,
        }
        url = BASE + "?" + urllib.parse.urlencode(params)
        try:
            req = urllib.request.Request(url, headers={"User-Agent":"feas-test"})
            with urllib.request.urlopen(req, timeout=60) as r:
                data = json.load(r)
        except Exception as e:
            print("  err", e, file=sys.stderr)
            time.sleep(3)
            pages += 1
            continue
        acts = data.get("activities", [])
        if not acts:
            break
        for a in acts:
            st = a.get("standard_type")
            if st not in ("IC50","Ki","Kd"):
                continue
            smi = a.get("canonical_smiles")
            mid = a.get("molecule_chembl_id")
            if smi and mid and mid not in mols:
                mols[mid] = smi
        pm = data.get("page_meta", {})
        tot = pm.get("total_count")
        offset += limit
        pages += 1
        print(f"  {tid}: page {pages}, uniq mols so far {len(mols)}, total_activities {tot}", file=sys.stderr)
        if not pm.get("next"):
            break
        time.sleep(0.5)
    return mols

all_data = {}
for name, tid in TARGETS.items():
    print("Fetching", name, tid, file=sys.stderr)
    mols = fetch_target(tid)
    all_data[name] = mols
    print(name, "->", len(mols), "unique active molecules", file=sys.stderr)

with open(os.path.join(OUT, "actives.json"), "w") as f:
    json.dump(all_data, f)
print("DONE")
for k,v in all_data.items():
    print(k, len(v))
