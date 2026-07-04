#!/usr/bin/env python
"""
BLIND HOLD-OUT scoring — R4 engine, identical backbone to score.py / exp_score.py.
Scores the externally-selected+labeled holdout_drugs.json against the CACHED 18-target
panel (ablation/abl_actives.json), z vs the SAME 24-drug background as score.py.
Strict LOO: remove query skeleton + skeletons of every other holdout drug sharing a
culprit target, from EVERY class. Recovery = best (min) rank across the drug's culprit_set.
"""
import json, os, statistics
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

HERE = os.path.dirname(os.path.abspath(__file__))
PANEL = json.load(open(os.path.join(HERE,"..","ablation","abl_actives.json")))  # cached ligands
HOLD  = json.load(open(os.path.join(HERE,"holdout_drugs.json")))

ASSAY = {
 "hERG_CHEMBL240":"hERG patch-clamp","Nav15_CHEMBL1980":"Nav1.5 electrophysiology",
 "Cav12_CHEMBL1940":"Cav1.2 electrophysiology","5HT2A_CHEMBL224":"5-HT2A binding",
 "5HT2B_CHEMBL1833":"5-HT2B counter-screen","5HT3_CHEMBL1899":"5-HT3 binding",
 "D2_CHEMBL217":"D2 binding","M2_CHEMBL211":"M2 muscarinic binding",
 "alpha1A_CHEMBL229":"alpha1A binding","H1_CHEMBL231":"H1 binding",
 "muOpioid_CHEMBL233":"mu-opioid binding","SERT_CHEMBL228":"SERT uptake",
 "NET_CHEMBL222":"NET uptake","DAT_CHEMBL238":"DAT uptake","CB1_CHEMBL218":"CB1 counter-screen",
 "AChE_CHEMBL220":"AChE inhibition","MAOA_CHEMBL1951":"MAO-A inhibition","COX1_CHEMBL221":"COX-1 inhibition",
}
# Baseline B: fixed, candidate-INDEPENDENT default order (same as exp_score.py / task spec)
BASELINE_B_ORDER = ["hERG_CHEMBL240","Nav15_CHEMBL1980","Cav12_CHEMBL1940","5HT2A_CHEMBL224",
 "D2_CHEMBL217","M2_CHEMBL211","alpha1A_CHEMBL229","H1_CHEMBL231","muOpioid_CHEMBL233",
 "SERT_CHEMBL228","5HT2B_CHEMBL1833","5HT3_CHEMBL1899","NET_CHEMBL222","DAT_CHEMBL238",
 "CB1_CHEMBL218","AChE_CHEMBL220","MAOA_CHEMBL1951","COX1_CHEMBL221"]

BACKGROUND = {
 "metformin":"CN(C)C(=N)N=C(N)N","aspirin":"CC(=O)OC1=CC=CC=C1C(=O)O",
 "acetaminophen":"CC(=O)NC1=CC=C(C=C1)O",
 "atorvastatin":"CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4",
 "omeprazole":"CC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=CC=C(C=C3N2)OC",
 "amoxicillin":"CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=C(C=C3)O)N)C(=O)O)C",
 "ibuprofen":"CC(C)CC1=CC=C(C=C1)C(C)C(=O)O","caffeine":"CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
 "lisinopril":"C1CC(N(C1)C(=O)C(CCCCN)NC(CCC2=CC=CC=C2)C(=O)O)C(=O)O",
 "warfarin":"CC(=O)CC(C1=CC=CC=C1)C2=C(C3=CC=CC=C3OC2=O)O",
 "furosemide":"C1=CC=C(C(=C1)COC2=CC(=C(C=C2)S(=O)(=O)N)Cl)C(=O)O",
 "hydrochlorothiazide":"C1=CC2=C(C=C1S(=O)(=O)N)S(=O)(=O)NCN2",
 "metoprolol":"CC(C)NCC(COC1=CC=C(C=C1)CCOC)O",
 "amlodipine":"CCOC(=O)C1=C(NC(=C(C1C2=CC=CC=C2Cl)C(=O)OC)C)COCCN",
 "simvastatin":"CCC(C)(C)C(=O)OC1CC(C=C2C1C(C(C=C2)C)CCC3CC(CC(=O)O3)O)C",
 "gabapentin":"C1CCC(CC1)(CC(=O)O)CN","sertraline":"CNC1CCC(C2=CC=CC=C12)C3=CC(=C(C=C3)Cl)Cl",
 "ciprofloxacin":"C1CC1N2C=C(C(=O)C3=CC(=C(C=C32)N4CCNCC4)F)C(=O)O",
 "prednisone":"CC12CC(=O)C3C(C1CCC2(C(=O)CO)O)CCC4=CC(=O)C=CC34C",
 "levothyroxine":"C1=C(C=C(C(=C1I)OC2=CC(=C(C(=C2)I)O)I)I)CC(C(=O)O)N",
 "diazepam":"CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3",
 "ranitidine":"CNC(=CNS(=O)(=O)C)NCCSCC1=CC=C(O1)CN(C)C",
 "loratadine":"CCOC(=O)N1CCC(=C2C3=C(CCC4=C2N=CC=C4)C=C(C=C3)Cl)CC1",
 "naproxen":"CC(C1=CC2=CC=C(C=C2C=C1)OC)C(=O)O",
 "clopidogrel":"COC(=O)C(C1=CC=CC=C1Cl)N2CCC3=C(C2)C=CS3",
}

def fp(smi):
    m=Chem.MolFromSmiles(smi)
    if m is None: return None,None
    return AllChem.GetMorganFingerprintAsBitVect(m,2,nBits=2048), Chem.MolToInchiKey(m)
def skel(ik): return ik.split("-")[0] if ik else None

# panel fps keyed by skeleton for LOO
TARGET_FPS={}
for t,mols in PANEL.items():
    lst=[]
    for mid,smi in mols.items():
        f,ik=fp(smi)
        if f is not None: lst.append((skel(ik),f))
    TARGET_FPS[t]=lst
targets=list(ASSAY.keys())

def class_score(qfp,t,excl_skel):
    sims=[DataStructs.TanimotoSimilarity(qfp,f) for sk,f in TARGET_FPS[t] if sk not in excl_skel]
    sims.sort(reverse=True)
    return statistics.mean(sims[:5]) if sims else 0.0

# background z (each bg drug excludes only its own skeleton)
BG={t:[] for t in targets}
for n,s in BACKGROUND.items():
    f,ik=fp(s)
    for t in targets: BG[t].append(class_score(f,t,{skel(ik)}))
bg_stats={t:(statistics.mean(BG[t]),statistics.pstdev(BG[t])) for t in targets}
def z(v,t):
    m,sd=bg_stats[t]; return (v-m)/sd if sd>0 else float('nan')

# holdout fps + LOO partner map (share any culprit target)
QF={}; QSK={}
for pid,d in HOLD.items():
    f,ik=fp(d["smiles"]); QF[pid]=f; QSK[pid]=skel(ik)
def loo_excl(pid):
    ex={QSK[pid]}
    cs=set(HOLD[pid]["culprit_set"])
    for o,od in HOLD.items():
        if o==pid: continue
        if cs & set(od["culprit_set"]): ex.add(QSK[o])
    return ex

def rank_of(order,key): return order.index(key)+1

results={"assay":ASSAY,"baseline_B_order":[ASSAY[t] for t in BASELINE_B_ORDER],
         "bg_stats":bg_stats,"per_drug":{}}
cats=["buried","cardiac"]
rec={c:{m:{"top1":0,"top3":0,"n":0} for m in ["ours","B"]} for c in cats}
ranksum={c:{m:0 for m in ["ours","B"]} for c in cats}
money={c:0 for c in cats}; money_drugs={c:[] for c in cats}

rows=[]
for pid,d in HOLD.items():
    ex=loo_excl(pid)
    zv={t:z(class_score(QF[pid],t,ex),t) for t in targets}
    ours_order=sorted(targets,key=lambda t:zv[t],reverse=True)
    cs=d["culprit_set"]; prim=d["primary_culprit"]; cat=d["category"]
    # recovery = best (min) rank across culprit_set
    r_ours=min(rank_of(ours_order,k) for k in cs)
    r_B   =min(rank_of(BASELINE_B_ORDER,k) for k in cs)
    r_ours_prim=rank_of(ours_order,prim); r_B_prim=rank_of(BASELINE_B_ORDER,prim)
    kz=max(zv[k] for k in cs)  # best culprit z
    results["per_drug"][pid]={"name":d["name"],"category":cat,"culprit_type":d["culprit_type"],
        "culprit_set":cs,"primary_culprit":prim,"warning_classes":d["warning_classes"],
        "rank_ours":r_ours,"rank_B":r_B,"rank_ours_primary":r_ours_prim,"rank_B_primary":r_B_prim,
        "killer_z":kz,"ours_top3":[ASSAY[t] for t in ours_order[:3]],
        "z_vector":{ASSAY[t]:round(zv[t],2) for t in targets}}
    rec[cat]["ours"]["n"]+=1; rec[cat]["B"]["n"]+=1
    rec[cat]["ours"]["top1"]+=r_ours==1; rec[cat]["ours"]["top3"]+=r_ours<=3
    rec[cat]["B"]["top1"]+=r_B==1; rec[cat]["B"]["top3"]+=r_B<=3
    ranksum[cat]["ours"]+=r_ours; ranksum[cat]["B"]+=r_B
    if r_ours<=3 and r_B>3: money[cat]+=1; money_drugs[cat].append(d["name"])
    rows.append((d["name"],cat,d["culprit_type"],prim.split("_")[0],",".join(k.split("_")[0] for k in cs),r_ours,r_B,kz))

meanrank={c:{m:(ranksum[c][m]/rec[c]["ours"]["n"] if rec[c]["ours"]["n"] else 0) for m in ["ours","B"]} for c in cats}
results["recovery_by_category"]=rec
results["mean_assays_to_culprit"]=meanrank
results["money_number"]=money; results["money_drugs"]=money_drugs

print(f"\n{'drug':22} {'cat':7} {'ctype':11} {'primary':9} {'culprit_set':16} {'OURS':>4}{'B':>4} {'z':>6}")
for r in sorted(rows,key=lambda x:(x[1],x[0])):
    print(f"{r[0][:22]:22} {r[1]:7} {r[2]:11} {r[3]:9} {r[4][:16]:16} {r[5]:>4}{r[6]:>4} {r[7]:>+6.2f}")

print("\n=== RECOVERY by category (best rank across culprit_set) ===")
for c in cats:
    n=rec[c]["ours"]["n"]
    if not n: continue
    print(f"[{c}] n={n}")
    print(f"  Top-1: Ours {rec[c]['ours']['top1']}/{n}  B {rec[c]['B']['top1']}/{n}")
    print(f"  Top-3: Ours {rec[c]['ours']['top3']}/{n}  B {rec[c]['B']['top3']}/{n}")
    print(f"  Mean assays-to-culprit: Ours {meanrank[c]['ours']:.2f}  B {meanrank[c]['B']:.2f}")
print("\n=== MONEY (Ours top3 AND B not top3) ===")
for c in cats: print(f"  [{c}] {money[c]}/{rec[c]['ours']['n']} -> {money_drugs[c]}")

# organ-anchored subset (culprit label externally tied to the withdrawal reason)
print("\n=== ORGAN-ANCHORED SUBSET (culprit_type=='organ') ===")
org=[pid for pid,d in HOLD.items() if d["culprit_type"]=="organ"]
for pid in org:
    p=results["per_drug"][pid]
    print(f"  {p['name'][:22]:22} {p['category']:7} set={','.join(k.split('_')[0] for k in p['culprit_set'])} "
          f"Ours={p['rank_ours']} B={p['rank_B']} z={p['killer_z']:+.2f}")

json.dump(results,open(os.path.join(HERE,"results.json"),"w"),indent=1)
print("\nwrote results.json")
