import json, urllib.request, urllib.parse, time, sys, os

OUT = os.path.dirname(os.path.abspath(__file__))
BASE = "https://www.ebi.ac.uk/chembl/api/data/activity.json"

# 18-target safety panel: name -> ChEMBL target id
TARGETS = {
    "hERG_CHEMBL240":     "CHEMBL240",    # hERG patch-clamp
    "Nav15_CHEMBL1980":   "CHEMBL1980",   # Nav1.5 electrophysiology
    "Cav12_CHEMBL1940":   "CHEMBL1940",   # Cav1.2 electrophysiology
    "5HT2A_CHEMBL224":    "CHEMBL224",    # 5-HT2A binding
    "5HT2B_CHEMBL1833":   "CHEMBL1833",   # 5-HT2B counter-screen
    "5HT3_CHEMBL1899":    "CHEMBL1899",   # 5-HT3 binding
    "D2_CHEMBL217":       "CHEMBL217",    # D2 binding
    "M2_CHEMBL211":       "CHEMBL211",    # M2 muscarinic binding
    "alpha1A_CHEMBL229":  "CHEMBL229",    # alpha1A binding
    "H1_CHEMBL231":       "CHEMBL231",    # H1 binding
    "muOpioid_CHEMBL233": "CHEMBL233",    # mu-opioid binding
    "SERT_CHEMBL228":     "CHEMBL228",    # SERT uptake
    "NET_CHEMBL222":      "CHEMBL222",    # NET uptake
    "DAT_CHEMBL238":      "CHEMBL238",    # DAT uptake
    "CB1_CHEMBL218":      "CHEMBL218",    # CB1 counter-screen
    "AChE_CHEMBL220":     "CHEMBL220",    # AChE inhibition
    "MAOA_CHEMBL1951":    "CHEMBL1951",   # MAO-A inhibition
    "COX1_CHEMBL221":     "CHEMBL221",    # COX-1 inhibition
}

def fetch_target(tid, want=900):
    mols = {}
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
            req = urllib.request.Request(url, headers={"User-Agent": "exp-test"})
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
            if st not in ("IC50", "Ki", "Kd"):
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

with open(os.path.join(OUT, "exp_actives.json"), "w") as f:
    json.dump(all_data, f)
print("DONE")
for k, v in all_data.items():
    print(k, len(v))
