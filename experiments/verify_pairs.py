"""Verify 'fingerprints miss it' demo pairs.
Fill COMPOUNDS with {name: SMILES}; PAIRS with [(a,b,shared_mechanism), ...].
Prints ECFP4 & MACCS Tanimoto (low = fingerprint search misses the link)
plus physchem descriptors (to show a property/pharmacophore channel connects them)."""
import sys
from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys, Descriptors, Crippen, rdMolDescriptors
from rdkit.Chem import DataStructs
from rdkit.Chem.MolStandardize import rdMolStandardize

COMPOUNDS = {}   # e.g. "terfenadine": "OC(...)..."
PAIRS = []       # e.g. ("terfenadine", "cisapride", "hERG block / QT")

def basic_n(mol):
    # crude cationic-amphiphile proxy: count basic (protonatable) nitrogens
    patt = Chem.MolFromSmarts("[$([NX3;H0,H1,H2;!$(NC=O);!$(N=*);!$([N+])]);!$(Nc)]")
    return len(mol.GetSubstructMatches(patt)) if patt else 0

def desc(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None, None
    d = {
        "MW": round(Descriptors.MolWt(m), 1),
        "cLogP": round(Crippen.MolLogP(m), 2),
        "TPSA": round(rdMolDescriptors.CalcTPSA(m), 1),
        "HBD": rdMolDescriptors.CalcNumHBD(m),
        "HBA": rdMolDescriptors.CalcNumHBA(m),
        "ArRings": rdMolDescriptors.CalcNumAromaticRings(m),
        "basicN": basic_n(m),
    }
    return m, d

def ecfp(m):
    return AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048)

def run():
    mols = {}
    print("=== Descriptors ===")
    for name, smi in COMPOUNDS.items():
        m, d = desc(smi)
        if m is None:
            print(f"  !! {name}: INVALID SMILES: {smi}")
            continue
        mols[name] = m
        print(f"  {name:22s} {d}")
    print("\n=== Pair similarity (ECFP4 r2 2048b Tanimoto | MACCS Tanimoto) ===")
    for a, b, mech in PAIRS:
        if a not in mols or b not in mols:
            print(f"  !! {a} / {b}: missing valid mol"); continue
        ea, eb = ecfp(mols[a]), ecfp(mols[b])
        t_ecfp = DataStructs.TanimotoSimilarity(ea, eb)
        ma, mb = MACCSkeys.GenMACCSKeys(mols[a]), MACCSkeys.GenMACCSKeys(mols[b])
        t_maccs = DataStructs.TanimotoSimilarity(ma, mb)
        verdict = "GOOD (fingerprints miss it)" if t_ecfp < 0.35 else ("marginal" if t_ecfp < 0.5 else "TOO SIMILAR (obvious analog)")
        print(f"  {a} <> {b}")
        print(f"     ECFP4={t_ecfp:.3f}  MACCS={t_maccs:.3f}  [{verdict}]  shared: {mech}")

if __name__ == "__main__":
    run()
