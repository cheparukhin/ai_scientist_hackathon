"""Fetch ChEMBL actives for the 4 panel targets missing from the sibling
assay_recovery cache, then merge with the sibling's 14 to build the full
18-target abl_actives.json. Same fetch pattern as fetch_actives.py /
ar_fetch.py (pChEMBL>=6, standard_type in IC50/Ki/Kd)."""
import json, urllib.request, urllib.parse, time, sys, os

HERE = os.path.dirname(os.path.abspath(__file__))
SIBLING = os.path.join(HERE, "..", "assay_recovery", "ar_actives.json")
BASE = "https://www.ebi.ac.uk/chembl/api/data/activity.json"

# 4 targets in the 18-panel that the sibling cache does NOT have
NEW_TARGETS = {
    "5HT3_CHEMBL1899": "CHEMBL1899",   # 5-HT3 counter-screen
    "NET_CHEMBL222":   "CHEMBL222",    # norepinephrine transporter
    "DAT_CHEMBL238":   "CHEMBL238",    # dopamine transporter
    "CB1_CHEMBL218":   "CHEMBL218",    # cannabinoid CB1
}

def fetch_target(tid, want=900):
    mols = {}
    offset = 0
    limit = 1000
    pages = 0
    while len(mols) < want and pages < 8:
        params = {"target_chembl_id": tid, "pchembl_value__gte": 6,
                  "limit": limit, "offset": offset}
        url = BASE + "?" + urllib.parse.urlencode(params)
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "abl-test"})
            with urllib.request.urlopen(req, timeout=60) as r:
                data = json.load(r)
        except Exception as e:
            print("  err", e, file=sys.stderr); time.sleep(3); pages += 1; continue
        acts = data.get("activities", [])
        if not acts:
            break
        for a in acts:
            if a.get("standard_type") not in ("IC50", "Ki", "Kd"):
                continue
            smi = a.get("canonical_smiles"); mid = a.get("molecule_chembl_id")
            if smi and mid and mid not in mols:
                mols[mid] = smi
        pm = data.get("page_meta", {})
        offset += limit; pages += 1
        print(f"  {tid}: page {pages}, uniq {len(mols)}, total {pm.get('total_count')}", file=sys.stderr)
        if not pm.get("next"):
            break
        time.sleep(0.5)
    return mols

all_data = json.load(open(SIBLING))   # start from sibling's 14 targets
print("Loaded sibling cache:", {k: len(v) for k, v in all_data.items()}, file=sys.stderr)
for name, tid in NEW_TARGETS.items():
    print("Fetching", name, tid, file=sys.stderr)
    all_data[name] = fetch_target(tid)
    print(name, "->", len(all_data[name]), file=sys.stderr)

with open(os.path.join(HERE, "abl_actives.json"), "w") as f:
    json.dump(all_data, f)
print("DONE. 18-target counts:")
for k, v in all_data.items():
    print(k, len(v))
