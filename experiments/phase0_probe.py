#!/usr/bin/env python3
"""Phase-0 automatability probe: SMILES + ChEMBL safety-target bioactivity."""
import urllib.request, urllib.parse, json, time, sys

DRUGS = ["terfenadine","cisapride","astemizole","thioridazine","sertindole",
"grepafloxacin","cerivastatin","pergolide","fenfluramine","benfluorex",
"troglitazone","nefazodone","sitaxentan","fialuridine","perhexiline",
"rofecoxib","ximelagatran","tolcapone","trovafloxacin","bromfenac"]

# Curated safety / off-target panel (ChEMBL target IDs)
SAFETY_TARGETS = {
 "CHEMBL240":  "hERG (KCNH2) cardiac",
 "CHEMBL1833": "5-HT2B (HTR2B) valvulopathy",
 "CHEMBL224":  "5-HT2A",
 "CHEMBL6020": "BSEP (ABCB11) liver",
 "CHEMBL230":  "COX-2 (PTGS2)",
 "CHEMBL221":  "COX-1 (PTGS1)",
 "CHEMBL217":  "D2 receptor (DRD2)",
 "CHEMBL222":  "NET (SLC6A2)",
 "CHEMBL228":  "SERT (SLC6A4)",
 "CHEMBL251":  "Adenosine A2A",
 "CHEMBL2039": "MAO-A",
 "CHEMBL1951": "MAO-B",
 "CHEMBL4302": "Nav1.5 (SCN5A) cardiac",
 "CHEMBL4523987": "COMT",
 "CHEMBL3242": "Endothelin ET-A",
 "CHEMBL252":  "Endothelin ET-B",
 "CHEMBL4128": "Thrombin (F2)",
}

def fetch(url, tries=3):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent":"phase0probe/1.0","Accept":"application/json"})
            with urllib.request.urlopen(req, timeout=45) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            if i==tries-1:
                return {"__error__": str(e)}
            time.sleep(2)

def pubchem_smiles(name):
    base="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
    for prop in ["CanonicalSMILES","ConnectivitySMILES","SMILES"]:
        url=base+urllib.parse.quote(name)+f"/property/{prop}/JSON"
        d=fetch(url)
        if d and "PropertyTable" in d:
            props=d["PropertyTable"]["Properties"][0]
            for k in ("CanonicalSMILES","ConnectivitySMILES","SMILES"):
                if k in props: return props[k]
    return None

def chembl_molid(name):
    # search by pref_name exact, then by synonym
    url="https://www.ebi.ac.uk/chembl/api/data/molecule.json?pref_name__iexact="+urllib.parse.quote(name.upper())
    d=fetch(url)
    if d and d.get("molecules"):
        return d["molecules"][0]["molecule_chembl_id"]
    # synonym search
    url="https://www.ebi.ac.uk/chembl/api/data/molecule/search.json?q="+urllib.parse.quote(name)
    d=fetch(url)
    if d and d.get("molecules"):
        return d["molecules"][0]["molecule_chembl_id"]
    return None

def chembl_safety_activities(molid):
    """Return dict target_id -> count of measured activities (IC50/Ki/Kd/pChEMBL)."""
    hits={}
    tid_filter=",".join(SAFETY_TARGETS.keys())
    # ChEMBL activity endpoint filter by molecule + target set
    url=("https://www.ebi.ac.uk/chembl/api/data/activity.json?molecule_chembl_id="
         +molid+"&target_chembl_id__in="+urllib.parse.quote(tid_filter)+"&limit=1000")
    d=fetch(url)
    if not d or "activities" not in d:
        return hits, d.get("__error__") if d else "no data"
    for a in d["activities"]:
        t=a.get("target_chembl_id")
        if t in SAFETY_TARGETS:
            hits.setdefault(t,0)
            hits[t]+=1
    return hits, None

results={}
for name in DRUGS:
    print(f"=== {name} ===", file=sys.stderr)
    smiles=pubchem_smiles(name)
    molid=chembl_molid(name)
    safety={}
    err=None
    if molid:
        safety,err=chembl_safety_activities(molid)
    results[name]={
        "smiles": smiles,
        "smiles_ok": bool(smiles),
        "chembl_id": molid,
        "safety_targets": {SAFETY_TARGETS[t]:c for t,c in safety.items()},
        "n_safety_targets": len(safety),
        "err": err,
    }
    print(f"  smiles={'Y' if smiles else 'N'} chembl={molid} safety_targets={len(safety)} {list(results[name]['safety_targets'].keys())}", file=sys.stderr)
    time.sleep(0.3)

out="/private/tmp/claude-501/-Users-cheparukhin-hackathon/dd15a6b7-3182-4dbe-be6a-7324a3b185bb/scratchpad/phase0_results.json"
with open(out,"w") as f:
    json.dump(results,f,indent=2)

# Summary
n_smiles=sum(1 for r in results.values() if r["smiles_ok"])
n_chembl=sum(1 for r in results.values() if r["chembl_id"])
n_safety=sum(1 for r in results.values() if r["n_safety_targets"]>0)
print("\n===== SUMMARY =====")
print(f"SMILES resolved (PubChem): {n_smiles}/20")
print(f"ChEMBL molecule ID found:  {n_chembl}/20")
print(f">=1 measured safety-target activity in ChEMBL: {n_safety}/20")
print("\nPer-drug safety-target counts:")
for name in DRUGS:
    r=results[name]
    print(f"  {name:16s} chembl={str(r['chembl_id']):16s} nTargets={r['n_safety_targets']:2d}  {', '.join(r['safety_targets'].keys())}")
