"""build_outcomes.py - validate/normalize the committed M2 curated data.

Asserts every SMILES in outcome_actives.json and every SMARTS in reactive_alerts.json
parses in RDKit, canonicalizes the SMILES in place (idempotent), and reports counts.

Run:  .venv/bin/python app/build_outcomes.py
"""
import json
import os
import sys

from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")


def validate(write=False):
    ok = True
    actives = json.load(open(os.path.join(DATA, "outcome_actives.json")))
    alerts = json.load(open(os.path.join(DATA, "reactive_alerts.json")))

    for endpoint, drugs in actives.items():
        for d in drugs:
            m = Chem.MolFromSmiles(d["smiles"])
            if m is None:
                ok = False
                print(f"  [FAIL smiles] {endpoint}: {d['name']}: {d['smiles']!r}")
            else:
                d["smiles"] = Chem.MolToSmiles(m)
        print(f"  {endpoint}: {len(drugs)} drugs")

    for a in alerts:
        patt = Chem.MolFromSmarts(a["smarts"])
        if patt is None:
            ok = False
            print(f"  [FAIL smarts] {a['name']}: {a['smarts']!r}")
    print(f"  reactive_alerts: {len(alerts)} SMARTS")

    if ok and write:
        json.dump(actives, open(os.path.join(DATA, "outcome_actives.json"), "w"), indent=1)
        print("  (canonicalized SMILES written back)")

    return ok


if __name__ == "__main__":
    print("Validating M2 curated data...")
    ok = validate(write="--write" in sys.argv)
    if ok:
        print("All SMILES/SMARTS parse. OK.")
    else:
        print("PARSE FAILURES above - fix before committing.")
        sys.exit(1)
