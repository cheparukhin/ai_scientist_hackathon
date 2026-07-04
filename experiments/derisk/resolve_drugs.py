"""Resolve canonical SMILES + InChIKey for the expanded worth-it query set via PubChem.
I (Claude) own the toxicology: drug -> culprit off-target -> organ. PubChem owns structure.
Writes ../derisk/drugs.json for the parallel scoring/ablation experiments to consume.
"""
import json, os, time, urllib.request, urllib.parse
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

HERE = os.path.dirname(os.path.abspath(__file__))

# name -> (culprit_target_key, culprit_chembl_id, organ/phenotype, category)
# category: "buried" = culprit is NOT a front-loaded cardiac ion channel (the hypothesis test)
#           "cardiac_herg" = hERG/QT (expected no-op vs a default panel; contrast set)
DRUGS = {
    # --- buried-target failures: culprit is a non-cardiac-channel off-target ---
    "pergolide":       ("5-HT2B", "CHEMBL1833", "valvular fibrosis", "buried"),
    "cabergoline":     ("5-HT2B", "CHEMBL1833", "valvular fibrosis", "buried"),
    "fenfluramine":    ("5-HT2B", "CHEMBL1833", "valvulopathy (metabolite-active)", "buried"),
    "dexfenfluramine": ("5-HT2B", "CHEMBL1833", "valvulopathy (metabolite-active)", "buried"),
    "methysergide":    ("5-HT2B", "CHEMBL1833", "retroperitoneal/valve fibrosis", "buried"),
    "rimonabant":      ("CB1",    "CHEMBL218",  "psychiatric / suicidality", "buried"),
    "taranabant":      ("CB1",    "CHEMBL218",  "psychiatric", "buried"),
    "alosetron":       ("5-HT3",  "CHEMBL1899", "ischemic colitis", "buried"),
    "sibutramine":     ("NET",    "CHEMBL222",  "cardiovascular (BP/HR)", "buried"),
    "mibefradil":      ("Cav1.2", "CHEMBL1940", "cardiac (T-type Ca) + DDI", "buried"),
    # --- hERG/QT cardiac contrast set: default panels already front-load hERG ---
    "terfenadine":     ("hERG", "CHEMBL240", "QT / TdP", "cardiac_herg"),
    "cisapride":       ("hERG", "CHEMBL240", "QT / TdP", "cardiac_herg"),
    "astemizole":      ("hERG", "CHEMBL240", "QT / TdP", "cardiac_herg"),
    "grepafloxacin":   ("hERG", "CHEMBL240", "QT / TdP", "cardiac_herg"),
    "sertindole":      ("hERG", "CHEMBL240", "QT / TdP", "cardiac_herg"),
    "thioridazine":    ("hERG", "CHEMBL240", "QT / TdP", "cardiac_herg"),
    "terodiline":      ("hERG", "CHEMBL240", "QT / TdP", "cardiac_herg"),
    "droperidol":      ("hERG", "CHEMBL240", "QT / TdP", "cardiac_herg"),
    "sparfloxacin":    ("hERG", "CHEMBL240", "QT / TdP", "cardiac_herg"),
    "mesoridazine":    ("hERG", "CHEMBL240", "QT / TdP", "cardiac_herg"),
}

def pubchem_smiles(name):
    base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
    url = base + urllib.parse.quote(name) + "/property/SMILES,ConnectivitySMILES/JSON"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "worth-it-derisk"})
        with urllib.request.urlopen(req, timeout=30) as r:
            d = json.load(r)
        props = d["PropertyTable"]["Properties"][0]
        return props.get("SMILES") or props.get("ConnectivitySMILES")
    except Exception as e:
        print(f"  {name}: fail {e}")
        return None

out = {}
for name, (tkey, tid, organ, cat) in DRUGS.items():
    smi = pubchem_smiles(name)
    ik = None
    if smi:
        m = Chem.MolFromSmiles(smi)
        if m:
            smi = Chem.MolToSmiles(m)  # canonicalize with RDKit for consistency
            ik = Chem.MolToInchiKey(m)
        else:
            print(f"  {name}: RDKit failed to parse {smi}")
            smi = None
    out[name] = {"smiles": smi, "inchikey": ik, "culprit_target": tkey,
                 "culprit_chembl_id": tid, "organ": organ, "category": cat}
    print(f"{'OK ' if smi else 'FAIL'} {name:16s} {tkey:8s} {cat:12s} {ik or ''}")
    time.sleep(0.3)

json.dump(out, open(os.path.join(HERE, "drugs.json"), "w"), indent=1)
n_ok = sum(1 for v in out.values() if v["smiles"])
print(f"\nResolved {n_ok}/{len(DRUGS)} -> drugs.json")
print("buried:", sum(1 for v in out.values() if v["smiles"] and v["category"]=="buried"),
      "| cardiac_herg:", sum(1 for v in out.values() if v["smiles"] and v["category"]=="cardiac_herg"))
