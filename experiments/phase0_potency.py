#!/usr/bin/env python3
"""Refine: capture best pChEMBL per safety target -> distinguish potent hits from panel noise."""
import urllib.request, urllib.parse, json, time, sys

RES="/private/tmp/claude-501/-Users-cheparukhin-hackathon/dd15a6b7-3182-4dbe-be6a-7324a3b185bb/scratchpad/phase0_results.json"
prev=json.load(open(RES))

SAFETY_TARGETS={"CHEMBL240":"hERG","CHEMBL1833":"5-HT2B","CHEMBL224":"5-HT2A",
"CHEMBL6020":"BSEP","CHEMBL230":"COX-2","CHEMBL221":"COX-1","CHEMBL217":"D2",
"CHEMBL222":"NET","CHEMBL228":"SERT","CHEMBL251":"A2A","CHEMBL2039":"MAO-A",
"CHEMBL1951":"MAO-B","CHEMBL4302":"Nav1.5","CHEMBL4523987":"COMT",
"CHEMBL3242":"ET-A","CHEMBL252":"ET-B","CHEMBL4128":"Thrombin"}

def fetch(url,tries=3):
    for i in range(tries):
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"p/1.0","Accept":"application/json"})
            with urllib.request.urlopen(req,timeout=45) as r: return json.loads(r.read().decode())
        except Exception as e:
            if i==tries-1: return {"__error__":str(e)}
            time.sleep(2)

def best(molid):
    """target -> (best_pchembl, standard_type, n_records)"""
    d={}
    tid=",".join(SAFETY_TARGETS)
    url=("https://www.ebi.ac.uk/chembl/api/data/activity.json?molecule_chembl_id="+molid
         +"&target_chembl_id__in="+urllib.parse.quote(tid)+"&limit=1000")
    r=fetch(url)
    if not r or "activities" not in r: return d
    for a in r["activities"]:
        t=a.get("target_chembl_id")
        if t not in SAFETY_TARGETS: continue
        name=SAFETY_TARGETS[t]
        p=a.get("pchembl_value")
        st=a.get("standard_type")
        try: pv=float(p) if p else None
        except: pv=None
        cur=d.get(name,{"pchembl":None,"type":st,"n":0})
        cur["n"]+=1
        if pv is not None and (cur["pchembl"] is None or pv>cur["pchembl"]):
            cur["pchembl"]=pv; cur["type"]=st
        d[name]=cur
    return d

summary={}
for name,rec in prev.items():
    mid=rec["chembl_id"]
    b=best(mid) if mid else {}
    # potent = pchembl>=6 (<=1uM) on a binding/functional assay
    potent=[k for k,v in b.items() if v["pchembl"] is not None and v["pchembl"]>=6.0]
    measured_any=[k for k,v in b.items() if v["pchembl"] is not None]  # has a real quantitative value
    summary[name]={"chembl":mid,"targets":b,
        "n_any_record":len(b),"n_quantified":len(measured_any),"n_potent_1uM":len(potent),
        "potent_list":sorted(potent)}
    top=sorted([(v["pchembl"],k) for k,v in b.items() if v["pchembl"] is not None],reverse=True)[:5]
    print(f"{name:16s} anyRec={len(b):2d} quantified={len(measured_any):2d} potent<=1uM={len(potent):2d}  top: "
          +", ".join(f"{k}={pv:.1f}" for pv,k in top),file=sys.stderr)
    time.sleep(0.3)

json.dump(summary,open("/private/tmp/claude-501/-Users-cheparukhin-hackathon/dd15a6b7-3182-4dbe-be6a-7324a3b185bb/scratchpad/phase0_potency.json","w"),indent=2)

n_any=sum(1 for s in summary.values() if s["n_any_record"]>0)
n_q=sum(1 for s in summary.values() if s["n_quantified"]>0)
n_p=sum(1 for s in summary.values() if s["n_potent_1uM"]>0)
print("\n===== POTENCY-FILTERED SUMMARY =====")
print(f">=1 safety-target RECORD (any, incl weak/panel):     {n_any}/20")
print(f">=1 safety-target with a QUANTIFIED pChEMBL value:   {n_q}/20")
print(f">=1 POTENT safety off-target (<=1uM, pChEMBL>=6):    {n_p}/20")
