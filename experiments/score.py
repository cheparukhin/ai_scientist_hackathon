import json, os, statistics
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

HERE=os.path.dirname(os.path.abspath(__file__))
data=json.load(open(os.path.join(HERE,"actives.json")))

QUERIES={
 "terfenadine":  "CC(C)(C)C1=CC=C(C=C1)C(CCCN2CCC(CC2)C(C3=CC=CC=C3)(C4=CC=CC=C4)O)O",
 "cisapride":    "CO[C@H]1CN(CC[C@H]1NC(=O)C2=CC(=C(C=C2OC)N)Cl)CCCOC3=CC=C(C=C3)F",
 "astemizole":   "COC1=CC=C(C=C1)CCN2CCC(CC2)NC3=NC4=CC=CC=C4N3CC5=CC=C(C=C5)F",
 "thioridazine": "CN1CCCCC1CCN2C3=CC=CC=C3SC4=C2C=C(C=C4)SC",
 "pergolide":    "CCCN1C[C@@H](C[C@H]2[C@H]1CC3=CNC4=CC=CC2=C34)CSC",
 "fenfluramine": "CCNC(C)CC1=CC(=CC=C1)C(F)(F)F",
}
TRUE={"terfenadine":"hERG_CHEMBL240","cisapride":"hERG_CHEMBL240",
      "astemizole":"hERG_CHEMBL240","thioridazine":"hERG_CHEMBL240",
      "pergolide":"5HT2B_CHEMBL1833","fenfluramine":"5HT2B_CHEMBL1833"}
# mechanistic partners to also remove during LOO (co-remove known hERG/5HT2B partners)
HERG_QUERIES={"terfenadine","cisapride","astemizole","thioridazine"}
HT2B_QUERIES={"pergolide","fenfluramine"}

BACKGROUND={
 "metformin":"CN(C)C(=N)N=C(N)N",
 "aspirin":"CC(=O)OC1=CC=CC=C1C(=O)O",
 "acetaminophen":"CC(=O)NC1=CC=C(C=C1)O",
 "atorvastatin":"CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4",
 "omeprazole":"CC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=CC=C(C=C3N2)OC",
 "amoxicillin":"CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=C(C=C3)O)N)C(=O)O)C",
 "ibuprofen":"CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
 "caffeine":"CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
 "lisinopril":"C1CC(N(C1)C(=O)C(CCCCN)NC(CCC2=CC=CC=C2)C(=O)O)C(=O)O",
 "warfarin":"CC(=O)CC(C1=CC=CC=C1)C2=C(C3=CC=CC=C3OC2=O)O",
 "furosemide":"C1=CC=C(C(=C1)COC2=CC(=C(C=C2)S(=O)(=O)N)Cl)C(=O)O",
 "hydrochlorothiazide":"C1=CC2=C(C=C1S(=O)(=O)N)S(=O)(=O)NCN2",
 "metoprolol":"CC(C)NCC(COC1=CC=C(C=C1)CCOC)O",
 "amlodipine":"CCOC(=O)C1=C(NC(=C(C1C2=CC=CC=C2Cl)C(=O)OC)C)COCCN",
 "simvastatin":"CCC(C)(C)C(=O)OC1CC(C=C2C1C(C(C=C2)C)CCC3CC(CC(=O)O3)O)C",
 "gabapentin":"C1CCC(CC1)(CC(=O)O)CN",
 "sertraline":"CNC1CCC(C2=CC=CC=C12)C3=CC(=C(C=C3)Cl)Cl",
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
    ik=Chem.MolToInchiKey(m)
    f=AllChem.GetMorganFingerprintAsBitVect(m,2,nBits=2048)
    return f,ik

# Precompute target active fps + inchikeys
print("Building target fingerprints...")
TARGET_FPS={}
for tname,mols in data.items():
    lst=[]
    for mid,smi in mols.items():
        f,ik=fp(smi)
        if f is not None:
            lst.append((ik,f))
    TARGET_FPS[tname]=lst
    print(f"  {tname}: {len(lst)} valid fps")

def class_score(qfp, tname, exclude_iks):
    sims=[]
    for ik,f in TARGET_FPS[tname]:
        if ik in exclude_iks: continue
        sims.append(DataStructs.TanimotoSimilarity(qfp,f))
    sims.sort(reverse=True)
    top5=sims[:5]
    return statistics.mean(top5), max(sims), len(sims)

# query inchikeys for LOO exclusion
QIK={n:fp(s)[1] for n,s in QUERIES.items()}

def loo_exclude(qname):
    # remove the query itself + its mechanistic partners in same class
    ex=set()
    ex.add(QIK[qname])
    grp = HERG_QUERIES if qname in HERG_QUERIES else HT2B_QUERIES
    for g in grp:
        ex.add(QIK[g])
    return ex

targets=list(data.keys())
print("\nTargets order:", targets)

# Background scores per target (no exclusion needed; but exclude bg self trivially not in set anyway)
BG={t:[] for t in targets}
for n,s in BACKGROUND.items():
    f,ik=fp(s)
    if f is None:
        print("bad bg",n); continue
    for t in targets:
        mean5,mx,ncomp=class_score(f,t,set([ik]))
        BG[t].append((n,mean5,mx))

bg_stats={}
for t in targets:
    vals=[x[1] for x in BG[t]]
    bg_stats[t]=(statistics.mean(vals), statistics.pstdev(vals), max(vals))

# Query scores
rows={}
for qn,qs in QUERIES.items():
    qf,_=fp(qs)
    ex=loo_exclude(qn)
    r={}
    for t in targets:
        mean5,mx,ncomp=class_score(qf,t,ex)
        r[t]=(mean5,mx,ncomp)
    rows[qn]=r

# pairwise terfenadine-cisapride
tf=fp(QUERIES["terfenadine"])[0]; cf=fp(QUERIES["cisapride"])[0]
pw=DataStructs.TanimotoSimilarity(tf,cf)

out={"targets":targets,"bg_stats":bg_stats,"rows":rows,"pairwise_terf_cis":pw,
     "bg_detail":BG,"n_actives":{t:len(TARGET_FPS[t]) for t in targets}}
json.dump(out,open(os.path.join(HERE,"results.json"),"w"),indent=1)

# Print report table
print("\n=== PAIRWISE BASELINE ===")
print(f"terfenadine-cisapride ECFP4 Tanimoto = {pw:.3f}")

print("\n=== BACKGROUND (24 unrelated drugs), top5-mean class-score ===")
for t in targets:
    m,sd,mx=bg_stats[t]
    print(f"  {t}: mean={m:.3f} sd={sd:.3f} max={mx:.3f}")

def zscore(v,t):
    m,sd,_=bg_stats[t]
    return (v-m)/sd if sd>0 else float('nan')

print("\n=== QUERY CLASS-SCORES (top5-mean [max]) with z vs background ===")
hdr="drug".ljust(13)
for t in targets: hdr+= t.split('_')[0].rjust(14)
hdr+="  true".rjust(8)+" top?".rjust(6)
print(hdr)
correct=0
for qn in QUERIES:
    line=qn.ljust(13)
    scores={t:rows[qn][t][0] for t in targets}
    best=max(scores,key=scores.get)
    tt=TRUE[qn]
    for t in targets:
        v=rows[qn][t][0]; z=zscore(v,t)
        tag="*" if t==tt else " "
        line+=f"{v:.2f}(z{z:+.1f}){tag}".rjust(14)
    istop = best==tt
    if istop: correct+=1
    line+=("  "+tt.split('_')[0]).rjust(8)+("  Y" if istop else "  N").rjust(6)
    print(line)
print(f"\nTrue-target top-ranked: {correct}/{len(QUERIES)}")
