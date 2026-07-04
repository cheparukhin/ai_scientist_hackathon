import itertools
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit import DataStructs
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

compounds = {
"pergolide":      "CCCN1C[C@@H](C[C@H]2[C@H]1CC3=CNC4=CC=CC2=C34)CSC",
"fenfluramine":   "CCNC(C)CC1=CC(=CC=C1)C(F)(F)F",
"norfenfluramine":"CC(CC1=CC(=CC=C1)C(F)(F)F)N",
"cabergoline":    "CCNC(=O)N(CCCN(C)C)C(=O)[C@@H]1C[C@H]2[C@@H](CC3=CNC4=CC=CC2=C34)N(C1)CC=C",
"benfluorex":     "CC(CC1=CC(=CC=C1)C(F)(F)F)NCCOC(=O)C2=CC=CC=C2",
"terfenadine":    "CC(C)(C)C1=CC=C(C=C1)C(CCCN2CCC(CC2)C(C3=CC=CC=C3)(C4=CC=CC=C4)O)O",
"cisapride":      "CO[C@H]1CN(CC[C@H]1NC(=O)C2=CC(=C(C=C2OC)N)Cl)CCCOC3=CC=C(C=C3)F",
"astemizole":     "COC1=CC=C(C=C1)CCN2CCC(CC2)NC3=NC4=CC=CC=C4N3CC5=CC=C(C=C5)F",
"thioridazine":   "CN1CCCCC1CCN2C3=CC=CC=C3SC4=C2C=C(C=C4)SC",
"sertindole":     "C1CN(CCC1C2=CN(C3=C2C=C(C=C3)Cl)C4=CC=C(C=C4)F)CCN5CCNC5=O",
"grepafloxacin":  "CC1CN(CCN1)C2=C(C(=C3C(=C2)N(C=C(C3=O)C(=O)O)C4CC4)C)F",
"troglitazone":   "CC1=C(C2=C(CCC(O2)(C)COC3=CC=C(C=C3)CC4C(=O)NC(=O)S4)C(=C1O)C)C",
"nefazodone":     "CCC1=NN(C(=O)N1CCOC2=CC=CC=C2)CCCN3CCN(CC3)C4=CC(=CC=C4)Cl",
"sitaxentan":     "CC1=CC2=C(C=C1CC(=O)C3=C(C=CS3)S(=O)(=O)NC4=C(C(=NO4)C)Cl)OCO2",
"fialuridine":    "C1=C(C(=O)NC(=O)N1[C@H]2[C@H]([C@@H]([C@H](O2)CO)O)F)I",
"perhexiline":    "C1CCC(CC1)C(CC2CCCCN2)C3CCCCC3",
}

demo_pairs = [
("terfenadine","cisapride","hERG"),
("astemizole","thioridazine","hERG"),
("pergolide","fenfluramine","5-HT2B"),
("cabergoline","benfluorex","5-HT2B"),
("troglitazone","nefazodone","BSEP+mito"),
("troglitazone","sitaxentan","BSEP"),
("fialuridine","perhexiline","mito"),
]
sanity_pairs = [
("pergolide","cabergoline"),
("fenfluramine","norfenfluramine"),
]

NCONF = 30
SEED = 42

mols = {}
usrcat = {}   # name -> list of per-conformer USRCAT descriptors
ecfp = {}
failures = []

for name, smi in compounds.items():
    m = Chem.MolFromSmiles(smi)
    if m is None:
        failures.append((name, "SMILES parse failed"))
        continue
    ecfp[name] = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048)
    mh = Chem.AddHs(m)
    params = AllChem.ETKDGv3()
    params.randomSeed = SEED
    cids = AllChem.EmbedMultipleConfs(mh, numConfs=NCONF, params=params)
    if len(cids) == 0:
        failures.append((name, "ETKDG embedding produced 0 conformers"))
        continue
    try:
        AllChem.MMFFOptimizeMoleculeConfs(mh)
    except Exception as e:
        failures.append((name, f"MMFF opt warning: {e}"))
    descs = []
    for cid in cids:
        try:
            descs.append(rdMolDescriptors.GetUSRCAT(mh, confId=cid))
        except Exception:
            pass
    if not descs:
        failures.append((name, "USRCAT failed for all conformers"))
        continue
    mols[name] = mh
    usrcat[name] = descs
    print(f"{name:16s} confs={len(cids):3d} usrcat_ok={len(descs)}")

print("\n=== Failures ===")
if failures:
    for n, r in failures:
        print(f"  {n}: {r}")
else:
    print("  none")

def best_usrcat(a, b):
    best = 0.0
    for da in usrcat[a]:
        for db in usrcat[b]:
            s = rdMolDescriptors.GetUSRScore(da, db)
            if s > best:
                best = s
    return best

def ecfp_tan(a, b):
    return DataStructs.TanimotoSimilarity(ecfp[a], ecfp[b])

names = [n for n in compounds if n in usrcat]

# Precompute full pairwise matrices for ranking
ecfp_mat = {a:{} for a in names}
usr_mat = {a:{} for a in names}
for a, b in itertools.combinations(names, 2):
    e = ecfp_tan(a, b)
    u = best_usrcat(a, b)
    ecfp_mat[a][b] = ecfp_mat[b][a] = e
    usr_mat[a][b] = usr_mat[b][a] = u

def rank_of_partner(query, partner, mat):
    # rank all others by similarity descending; return rank of partner (1=most similar)
    others = [(n, mat[query][n]) for n in names if n != query]
    others.sort(key=lambda x: -x[1])
    for i, (n, s) in enumerate(others, 1):
        if n == partner:
            return i, s
    return None, None

print("\n=== DEMO PAIRS ===")
print(f"{'pair':32s} {'mech':10s} {'ECFP4':>7s} {'USRCAT':>7s} {'E-rk(a->b/b->a)':>16s} {'U-rk(a->b/b->a)':>16s}")
demo_rows = []
for a, b, mech in demo_pairs:
    e = ecfp_mat[a][b]
    u = usr_mat[a][b]
    er_ab,_ = rank_of_partner(a, b, ecfp_mat)
    er_ba,_ = rank_of_partner(b, a, ecfp_mat)
    ur_ab,_ = rank_of_partner(a, b, usr_mat)
    ur_ba,_ = rank_of_partner(b, a, usr_mat)
    demo_rows.append((f"{a}-{b}", mech, e, u, er_ab, er_ba, ur_ab, ur_ba))
    print(f"{a+'-'+b:32s} {mech:10s} {e:7.3f} {u:7.3f} {str(er_ab)+'/'+str(er_ba):>16s} {str(ur_ab)+'/'+str(ur_ba):>16s}")

print("\n=== SANITY PAIRS (true analogs) ===")
print(f"{'pair':32s} {'ECFP4':>7s} {'USRCAT':>7s} {'E-rk(a->b/b->a)':>16s} {'U-rk(a->b/b->a)':>16s}")
for a, b in sanity_pairs:
    e = ecfp_mat[a][b]
    u = usr_mat[a][b]
    er_ab,_ = rank_of_partner(a, b, ecfp_mat)
    er_ba,_ = rank_of_partner(b, a, ecfp_mat)
    ur_ab,_ = rank_of_partner(a, b, usr_mat)
    ur_ba,_ = rank_of_partner(b, a, usr_mat)
    print(f"{a+'-'+b:32s} {e:7.3f} {u:7.3f} {str(er_ab)+'/'+str(er_ba):>16s} {str(ur_ab)+'/'+str(ur_ba):>16s}")

# Summary stats on ranking improvement (use best/min rank direction per pair)
print("\n=== RANK SUMMARY (best direction per demo pair) ===")
improved=0; total=0; top3_usr=0; top3_ecfp=0
for a,b,mech in demo_pairs:
    er = min(rank_of_partner(a,b,ecfp_mat)[0], rank_of_partner(b,a,ecfp_mat)[0])
    ur = min(rank_of_partner(a,b,usr_mat)[0], rank_of_partner(b,a,usr_mat)[0])
    total+=1
    if ur < er: improved+=1
    if ur<=3: top3_usr+=1
    if er<=3: top3_ecfp+=1
    print(f"{a}-{b:20s} ECFP4-rank={er}  USRCAT-rank={ur}  {'IMPROVED' if ur<er else ('WORSE' if ur>er else 'same')}")
print(f"\nDemo pairs where USRCAT rank < ECFP4 rank: {improved}/{total}")
print(f"Demo partners in top-3 neighbors: ECFP4={top3_ecfp}/{total}  USRCAT={top3_usr}/{total}")
print(f"(library size = {len(names)} compounds, so worst rank = {len(names)-1})")
