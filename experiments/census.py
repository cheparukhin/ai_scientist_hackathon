#!/usr/bin/env python
"""Data-coverage census for the Bowes 2012 secondary-pharmacology safety panel.
Counts distinct molecules with pChEMBL>=6 (IC50/Ki/Kd) per target in ChEMBL."""
import urllib.request, urllib.parse, json, time, sys

BASE = "https://www.ebi.ac.uk/chembl/api/data"
UA = {"User-Agent": "coverage-census/1.0"}

def get(url, tries=5):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.load(r)
        except Exception as e:
            sys.stderr.write(f"  retry {i} {e}\n")
            time.sleep(2*(i+1))
    raise RuntimeError("failed: "+url)

# target definitions: label -> (search query, preferred_chembl_id or None, organ/tox relevance)
# preferred id used to disambiguate; must be a SINGLE PROTEIN human target.
TARGETS = [
 # Ion channels
 ("hERG/KCNH2","hERG","CHEMBL240","Cardiac: QT prolongation, torsades"),
 ("Nav1.5/SCN5A","SCN5A","CHEMBL4296","Cardiac: conduction, arrhythmia"),
 ("Cav1.2/CACNA1C","CACNA1C",None,"Cardiac: contractility, BP"),
 ("KCNQ1","KCNQ1",None,"Cardiac: QT / LQT1"),
 # Aminergic GPCRs
 ("5-HT1A","serotonin 1a","CHEMBL214","CNS: mood, anxiety"),
 ("5-HT1B","serotonin 1b","CHEMBL1898","CNS/vascular"),
 ("5-HT2A","serotonin 2a","CHEMBL224","CNS: hallucination, antipsychotic"),
 ("5-HT2B","serotonin 2b","CHEMBL1833","Cardiac: valvulopathy, fibrosis"),
 ("5-HT3","serotonin 3a","CHEMBL339","GI: nausea, motility"),
 ("D1","dopamine d1","CHEMBL2056","CNS: motor, reward"),
 ("D2","dopamine d2","CHEMBL217","CNS: EPS, prolactin"),
 ("alpha1A adrenergic","adrenergic alpha-1a","CHEMBL229","Cardiovascular: BP, orthostasis"),
 ("alpha2A adrenergic","adrenergic alpha-2a","CHEMBL1867","CV/CNS: BP, sedation"),
 ("beta1 adrenergic","adrenergic beta-1","CHEMBL213","Cardiac: HR, contractility"),
 ("beta2 adrenergic","adrenergic beta-2","CHEMBL210","Respiratory/CV"),
 ("M1 muscarinic","muscarinic acetylcholine receptor m1","CHEMBL216","CNS/autonomic"),
 ("M2 muscarinic","muscarinic acetylcholine receptor m2","CHEMBL211","Cardiac: HR"),
 ("M3 muscarinic","muscarinic acetylcholine receptor m3","CHEMBL245","Smooth muscle, secretion"),
 ("H1 histamine","histamine h1","CHEMBL231","CNS: sedation"),
 ("H2 histamine","histamine h2","CHEMBL1250","GI: acid secretion"),
 ("mu opioid","mu opioid","CHEMBL233","CNS: analgesia, resp depression"),
 ("kappa opioid","kappa opioid","CHEMBL237","CNS: dysphoria"),
 ("delta opioid","delta opioid","CHEMBL236","CNS: analgesia, seizure"),
 ("CB1 cannabinoid","cannabinoid cb1","CHEMBL218","CNS: psychiatric"),
 ("CB2 cannabinoid","cannabinoid cb2","CHEMBL253","Immune"),
 ("A2A adenosine","adenosine a2a","CHEMBL251","CV/CNS"),
 # Transporters
 ("SERT/SLC6A4","serotonin transporter","CHEMBL228","CNS: serotonin syndrome"),
 ("NET/SLC6A2","norepinephrine transporter","CHEMBL222","CV/CNS: BP"),
 ("DAT/SLC6A3","dopamine transporter","CHEMBL238","CNS: abuse liability"),
 ("BSEP/ABCB11","bile salt export pump","CHEMBL6020","Liver: cholestatic DILI"),
 # Enzymes
 ("COX1/PTGS1","cyclooxygenase-1","CHEMBL221","GI: bleeding, ulcer"),
 ("COX2/PTGS2","cyclooxygenase-2","CHEMBL230","CV/renal"),
 ("AChE","acetylcholinesterase","CHEMBL220","Autonomic/CNS: cholinergic crisis"),
 ("MAO-A","monoamine oxidase a","CHEMBL1951","CNS/CV: hypertensive crisis"),
 ("PDE3A","phosphodiesterase 3a","CHEMBL4581","Cardiac: inotropy, mortality"),
 ("PDE4D","phosphodiesterase 4d","CHEMBL2622","CNS/GI: emesis"),
 ("Lck","lck","CHEMBL258","Immune: immunosuppression"),
 # Nuclear/other
 ("GR/NR3C1","glucocorticoid receptor","CHEMBL2034","Endocrine/metabolic"),
 ("androgen receptor","androgen receptor","CHEMBL1871","Endocrine/reproductive"),
]

ACT_TYPES = {"IC50","Ki","Kd"}

def resolve(query, pref):
    d = get(f"{BASE}/target/search.json?q={urllib.parse.quote(query)}&limit=20")
    cands = d.get("targets", [])
    if pref:
        for t in cands:
            if t["target_chembl_id"] == pref:
                return t
    # else pick best: human single protein
    human_sp = [t for t in cands if t.get("organism")=="Homo sapiens" and t.get("target_type")=="SINGLE PROTEIN"]
    if human_sp:
        # rank by score desc
        human_sp.sort(key=lambda t: -float(t.get("score",0)))
        return human_sp[0]
    if cands:
        return cands[0]
    return None

def count_actives(tid):
    """Distinct molecule_chembl_ids with pchembl>=6 and standard_type in IC50/Ki/Kd."""
    mols = set()
    total_records = 0
    url = (f"{BASE}/activity.json?target_chembl_id={tid}"
           f"&pchembl_value__gte=6&limit=1000"
           f"&only=molecule_chembl_id,standard_type,pchembl_value")
    while url:
        d = get(url)
        for a in d.get("activities", []):
            total_records += 1
            st = a.get("standard_type")
            if st in ACT_TYPES and a.get("molecule_chembl_id"):
                mols.add(a["molecule_chembl_id"])
        nxt = d.get("page_meta",{}).get("next")
        url = ("https://www.ebi.ac.uk"+nxt) if nxt else None
    return len(mols), total_records

def classify(n):
    if n >= 300: return "RICH"
    if n >= 50: return "USABLE"
    return "THIN"

results = []
for label, q, pref, relevance in TARGETS:
    sys.stderr.write(f"[{label}] resolving...\n")
    t = resolve(q, pref)
    if not t:
        results.append({"target":label,"chembl_id":None,"pref_name":None,"n_actives":None,"class":"UNRESOLVED","relevance":relevance})
        sys.stderr.write(f"  UNRESOLVED\n")
        continue
    tid = t["target_chembl_id"]
    n, recs = count_actives(tid)
    cls = classify(n)
    results.append({"target":label,"chembl_id":tid,"pref_name":t.get("pref_name"),
                    "organism":t.get("organism"),"target_type":t.get("target_type"),
                    "n_actives":n,"n_activity_records":recs,"class":cls,"relevance":relevance})
    sys.stderr.write(f"  {tid} {t.get('pref_name')} -> {n} actives ({recs} recs) [{cls}]\n")

with open("census_results.json","w") as f:
    json.dump(results, f, indent=2)

# CSV
import csv
with open("census_results.csv","w",newline="") as f:
    w = csv.writer(f)
    w.writerow(["target","chembl_id","pref_name","n_actives_pchembl6","class","organ_tox_relevance"])
    for r in results:
        w.writerow([r["target"],r["chembl_id"],r.get("pref_name"),r["n_actives"],r["class"],r["relevance"]])

print("DONE", len(results), "targets")
