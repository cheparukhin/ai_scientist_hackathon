from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit.Chem.Pharm2D import Generate, Gobbi_Pharm2D
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

C = {
  "pergolide":"CCCN1C[C@@H](C[C@H]2[C@H]1CC3=CNC4=CC=CC2=C34)CSC",
  "fenfluramine":"CCNC(C)CC1=CC(=CC=C1)C(F)(F)F",
  "cabergoline":"CCNC(=O)N(CCCN(C)C)C(=O)[C@@H]1C[C@H]2[C@@H](CC3=CNC4=CC=CC2=C34)N(C1)CC=C",
  "benfluorex":"CC(CC1=CC(=CC=C1)C(F)(F)F)NCCOC(=O)C2=CC=CC=C2",
  "terfenadine":"CC(C)(C)C1=CC=C(C=C1)C(CCCN2CCC(CC2)C(C3=CC=CC=C3)(C4=CC=CC=C4)O)O",
  "cisapride":"CO[C@H]1CN(CC[C@H]1NC(=O)C2=CC(=C(C=C2OC)N)Cl)CCCOC3=CC=C(C=C3)F",
  "astemizole":"COC1=CC=C(C=C1)CCN2CCC(CC2)NC3=NC4=CC=CC=C4N3CC5=CC=C(C=C5)F",
  "thioridazine":"CN1CCCCC1CCN2C3=CC=CC=C3SC4=C2C=C(C=C4)SC",
}
m = {k:Chem.MolFromSmiles(s) for k,s in C.items()}
fact = Gobbi_Pharm2D.factory
def ecfp(x): return AllChem.GetMorganFingerprintAsBitVect(x,2,nBits=2048)
def fcfp(x): return AllChem.GetMorganFingerprintAsBitVect(x,2,nBits=2048,useFeatures=True)
def gobbi(x): return Generate.Gen2DFingerprint(x, fact)
def t(f,a,b): return DataStructs.TanimotoSimilarity(f(m[a]),f(m[b]))

pairs = [
 ("pergolide","fenfluramine","5-HT2B valvulopathy"),
 ("cabergoline","benfluorex","5-HT2B valvulopathy"),
 ("terfenadine","cisapride","hERG/QT"),
 ("astemizole","thioridazine","hERG/QT"),
 # sanity: true analogs should stay high on all channels
 ("pergolide","cabergoline","SANITY ergolines"),
]
print(f"{'pair':32s} {'ECFP4':>7} {'FCFP4':>7} {'Gobbi':>7}   lift(best/ECFP)")
for a,b,lab in pairs:
    e,f,g = t(ecfp,a,b), t(fcfp,a,b), t(gobbi,a,b)
    best = max(f,g)
    lift = best/e if e>0 else float('inf')
    print(f"{a[:14]+'x'+b[:14]:32s} {e:7.3f} {f:7.3f} {g:7.3f}   {lift:4.1f}x   [{lab}]")
