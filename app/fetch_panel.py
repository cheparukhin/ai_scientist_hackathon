"""Build app/data/panel_actives.json from ChEMBL (18-target safety panel).

Reuses the exp_fetch.py logic: for each panel target, page ChEMBL activities at
pChEMBL >= 6, keep IC50/Ki/Kd, collect unique {molecule_chembl_id: canonical_smiles}.
Targets + ChEMBL IDs are read from app/data/panel.json.

Output panel_actives.json is a regenerable cache (gitignored). ~600-1500 actives/target.

Run:  .venv/bin/python app/fetch_panel.py
"""
import json
import os
import sys
import time
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = "https://www.ebi.ac.uk/chembl/api/data/activity.json"

PANEL = json.load(open(os.path.join(HERE, "data", "panel.json")))
TARGETS = {key: meta["chembl_id"] for key, meta in PANEL.items()}


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
            req = urllib.request.Request(url, headers={"User-Agent": "panel-fetch"})
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
            if a.get("standard_type") not in ("IC50", "Ki", "Kd"):
                continue
            smi = a.get("canonical_smiles")
            mid = a.get("molecule_chembl_id")
            if smi and mid and mid not in mols:
                mols[mid] = smi
        pm = data.get("page_meta", {})
        offset += limit
        pages += 1
        print(f"  {tid}: page {pages}, uniq mols so far {len(mols)}", file=sys.stderr)
        if not pm.get("next"):
            break
        time.sleep(0.5)
    return mols


def main():
    all_data = {}
    for name, tid in TARGETS.items():
        print("Fetching", name, tid, file=sys.stderr)
        mols = fetch_target(tid)
        all_data[name] = mols
        print(name, "->", len(mols), "unique active molecules", file=sys.stderr)
    out = os.path.join(HERE, "data", "panel_actives.json")
    with open(out, "w") as f:
        json.dump(all_data, f)
    print("DONE ->", out)
    for k, v in all_data.items():
        print(k, len(v))


if __name__ == "__main__":
    main()
