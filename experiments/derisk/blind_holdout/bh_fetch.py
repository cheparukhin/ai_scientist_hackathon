#!/usr/bin/env python
"""
BLIND HOLD-OUT fetch + externally-anchored labeling.

Every selection/labeling step is anchored to ChEMBL, NOT to us:
  1. SELECT drugs   = ChEMBL withdrawn drug_warnings (safety only).
  2. GET SMILES     = ChEMBL molecule.json parent structure.
  3. LABEL culprit  = ChEMBL measured activities (pChEMBL>=6, IC50/Ki/Kd) vs the 18-target panel.
Outputs holdout_drugs.json (externally selected + labeled, disjoint from our hand-picked drugs.json).
Run with: /Users/cheparukhin/hackathon/.venv/bin/python bh_fetch.py
"""
import urllib.request, urllib.parse, json, os, time, collections
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = "https://www.ebi.ac.uk/chembl/api/data/"

# ---- 18-target panel (matches cached abl_actives.json keys) ----
PANEL = {  # panel_key : chembl_target_id
    "hERG_CHEMBL240":"CHEMBL240","Nav15_CHEMBL1980":"CHEMBL1980","Cav12_CHEMBL1940":"CHEMBL1940",
    "5HT2A_CHEMBL224":"CHEMBL224","5HT2B_CHEMBL1833":"CHEMBL1833","5HT3_CHEMBL1899":"CHEMBL1899",
    "D2_CHEMBL217":"CHEMBL217","M2_CHEMBL211":"CHEMBL211","alpha1A_CHEMBL229":"CHEMBL229",
    "H1_CHEMBL231":"CHEMBL231","muOpioid_CHEMBL233":"CHEMBL233","SERT_CHEMBL228":"CHEMBL228",
    "NET_CHEMBL222":"CHEMBL222","DAT_CHEMBL238":"CHEMBL238","CB1_CHEMBL218":"CHEMBL218",
    "AChE_CHEMBL220":"CHEMBL220","MAOA_CHEMBL1951":"CHEMBL1951","COX1_CHEMBL221":"CHEMBL221",
}
TID2KEY = {v:k for k,v in PANEL.items()}
PANEL_IDS = ",".join(PANEL.values())
CARDIAC_CHANNELS = {"hERG_CHEMBL240","Nav15_CHEMBL1980","Cav12_CHEMBL1940"}

# ---- organ tags per panel target (from the panel doc) used for organ-matched culprit ----
def organ_targets(warning_classes, desc):
    """Panel targets whose organ matches the ChEMBL warning_class(es)/description."""
    desc = (desc or "").lower()
    S = set()
    for wc in warning_classes:
        if wc == "cardiotoxicity":
            if any(k in desc for k in ["valv","fibros","regurg","valve"]):
                S |= {"5HT2B_CHEMBL1833"}
            else:
                S |= {"hERG_CHEMBL240","Nav15_CHEMBL1980","Cav12_CHEMBL1940"}
        elif wc == "vascular toxicity":
            S |= {"alpha1A_CHEMBL229","NET_CHEMBL222"}
        elif wc == "psychiatric toxicity":
            S |= {"D2_CHEMBL217","CB1_CHEMBL218","SERT_CHEMBL228","5HT2A_CHEMBL224",
                  "DAT_CHEMBL238","MAOA_CHEMBL1951","NET_CHEMBL222"}
        elif wc == "neurotoxicity":
            S |= {"5HT2A_CHEMBL224","D2_CHEMBL217","H1_CHEMBL231","muOpioid_CHEMBL233",
                  "SERT_CHEMBL228","NET_CHEMBL222","DAT_CHEMBL238","CB1_CHEMBL218",
                  "AChE_CHEMBL220","MAOA_CHEMBL1951","M2_CHEMBL211"}
        elif wc == "gastrointestinal toxicity":
            S |= {"5HT3_CHEMBL1899","COX1_CHEMBL221","M2_CHEMBL211"}
        elif wc == "respiratory toxicity":
            S |= {"muOpioid_CHEMBL233","H1_CHEMBL231"}
        # hepatotoxicity / dermatological / immune / carcinogenicity / haematological /
        # nephrotoxicity / musculoskeletal / ocular / metabolic / terato / oto -> no panel organ
    return S

DROP_CLASSES = {None, "drug misuse"}

def get(url, tries=6):
    for t in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=90) as r:
                return json.load(r)
        except Exception as e:
            if t == tries-1: raise
            time.sleep(2*(t+1))

def skeleton(ik):
    return ik.split("-")[0] if ik else None

# ---------- STEP 1: withdrawn warnings (cached) ----------
raw_path = os.path.join(HERE,"raw_withdrawn.json")
if os.path.exists(raw_path):
    warnings = json.load(open(raw_path))
else:
    warnings=[]; url=BASE+"drug_warning.json?warning_type=Withdrawn&limit=1000"
    while url:
        d=get(url); warnings+=d["drug_warnings"]
        nxt=d["page_meta"].get("next"); url=("https://www.ebi.ac.uk"+nxt) if nxt else None
    json.dump(warnings, open(raw_path,"w"))
print(f"[1] withdrawn warning rows: {len(warnings)}")

# aggregate per parent
parents={}
for w in warnings:
    pid=w.get("parent_molecule_chembl_id")
    if not pid: continue
    p=parents.setdefault(pid,{"classes":set(),"descs":[]})
    p["classes"].add(w["warning_class"])
    if w.get("warning_description"): p["descs"].append(w["warning_description"])
print(f"[1] unique parents (all withdrawals): {len(parents)}")

# ---------- STEP 1b: SAFETY filter ----------
safety={}
for pid,p in parents.items():
    safety_classes={c for c in p["classes"] if c not in DROP_CLASSES}
    if not safety_classes: continue      # only None/misuse -> drop (efficacy/commercial/misuse)
    safety[pid]={"classes":sorted(safety_classes),"desc":" || ".join(p["descs"])}
print(f"[1b] parents with >=1 SAFETY withdrawal: {len(safety)}")

# ---------- STEP 3: measured activities vs panel (cached) ----------
act_path=os.path.join(HERE,"raw_activities.json")
if os.path.exists(act_path):
    ACT=json.load(open(act_path))
else:
    ACT={}
    ids=sorted(safety)
    for i,pid in enumerate(ids):
        url=(BASE+"activity.json?molecule_chembl_id=%s&target_chembl_id__in=%s"
             "&pchembl_value__gte=6&limit=1000"%(pid,PANEL_IDS))
        recs=[]
        while url:
            d=get(url); recs+=d["activities"]
            nxt=d["page_meta"].get("next"); url=("https://www.ebi.ac.uk"+nxt) if nxt else None
        ACT[pid]=recs
        if (i+1)%25==0: print(f"    activities {i+1}/{len(ids)}")
    json.dump(ACT, open(act_path,"w"))
print(f"[3] activity records fetched for {len(ACT)} parents")

# measured culprit set per parent: panel target with pChEMBL>=6 & type in IC50/Ki/Kd
OKTYPE={"IC50","Ki","Kd"}
measured={}  # pid -> {panel_key: max_pchembl}
for pid,recs in ACT.items():
    mp={}
    for r in recs:
        if r.get("standard_type") not in OKTYPE: continue
        pv=r.get("pchembl_value")
        if pv is None: continue
        try: pv=float(pv)
        except: continue
        if pv<6: continue
        tid=r.get("target_chembl_id")
        key=TID2KEY.get(tid)
        if not key: continue
        mp[key]=max(mp.get(key,0.0),pv)
    if mp: measured[pid]=mp
print(f"[4] parents with >=1 panel culprit (pChEMBL>=6): {len(measured)}")

# ---------- STEP 2: SMILES for surviving parents (cached) ----------
smi_path=os.path.join(HERE,"raw_smiles.json")
if os.path.exists(smi_path):
    SMI=json.load(open(smi_path))
else:
    SMI={}
    ids=sorted(measured)
    for i,pid in enumerate(ids):
        try:
            d=get(BASE+"molecule/%s.json"%pid)
            ms=d.get("molecule_structures") or {}
            SMI[pid]=ms.get("canonical_smiles")
            name=d.get("pref_name")
            SMI[pid+"__name"]=name
        except Exception:
            SMI[pid]=None
        if (i+1)%25==0: print(f"    smiles {i+1}/{len(ids)}")
    json.dump(SMI, open(smi_path,"w"))

# ---------- STEP 4 (disjoint) + build holdout records ----------
own=json.load(open(os.path.join(HERE,"..","drugs.json")))
OWN_SKEL={skeleton(v["inchikey"]) for v in own.values()}
print(f"[4] hand-picked drugs.json skeletons: {len(OWN_SKEL)}")

holdout={}; dropped_overlap=[]; dropped_nosmiles=[]
for pid,mp in measured.items():
    smi=SMI.get(pid)
    if not smi:
        dropped_nosmiles.append(pid); continue
    m=Chem.MolFromSmiles(smi)
    if m is None:
        dropped_nosmiles.append(pid); continue
    csmi=Chem.MolToSmiles(m)
    ik=Chem.MolToInchiKey(m)
    skel=skeleton(ik)
    if skel in OWN_SKEL:
        dropped_overlap.append((pid,SMI.get(pid+"__name"))); continue
    s=safety[pid]
    org=organ_targets(s["classes"], s["desc"])
    organ_culprits=sorted(set(mp) & org)
    if organ_culprits:
        culprit_set=organ_culprits; culprit_type="organ"
    else:
        best=max(mp, key=mp.get); culprit_set=[best]; culprit_type="max_pchembl"
    # single primary culprit = highest pChEMBL within culprit_set
    primary=max(culprit_set, key=lambda k: mp[k])
    category="cardiac" if all(k in CARDIAC_CHANNELS for k in culprit_set) else "buried"
    holdout[pid]={
        "chembl_id":pid,
        "name":SMI.get(pid+"__name"),
        "smiles":csmi,
        "inchikey":ik,
        "skeleton":skel,
        "warning_classes":s["classes"],
        "warning_desc":s["desc"][:400],
        "measured_culprits":{k:round(mp[k],2) for k in sorted(mp,key=lambda x:-mp[x])},
        "organ_target_set":sorted(org),
        "culprit_type":culprit_type,          # organ | max_pchembl
        "culprit_set":culprit_set,            # recovery = best rank across these
        "primary_culprit":primary,
        "max_pchembl_culprit":max(mp,key=mp.get),
        "category":category,
    }

print(f"[4] dropped (overlap w/ drugs.json): {len(dropped_overlap)} -> {dropped_overlap}")
print(f"[4] dropped (no/invalid SMILES): {len(dropped_nosmiles)}")
print(f"[FINAL] holdout drugs: {len(holdout)}")
cat=collections.Counter(v["category"] for v in holdout.values())
print(f"        by category: {dict(cat)}")
ctype=collections.Counter(v["culprit_type"] for v in holdout.values())
print(f"        culprit_type: {dict(ctype)}")

json.dump(holdout, open(os.path.join(HERE,"holdout_drugs.json"),"w"), indent=1)
print("wrote holdout_drugs.json")
