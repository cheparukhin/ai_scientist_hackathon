"""Probe which physicochemical descriptors separate OOD from in-domain drugs.
Computes descriptor ranges over the 7093-active reference set (OECD descriptor
domain) and compares in-domain drugs vs OOD on each axis, to build a defensible
property-bounding-box AD rule (since NN-distance alone failed to separate)."""
import json, os, statistics
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors, Crippen
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
HERE = os.path.dirname(os.path.abspath(__file__))
actives = json.load(open(os.path.join(HERE, "ad_actives.json")))

from calibrate_ad import QUERIES, BACKGROUND, OOD

def desc(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    heavy = m.GetNumHeavyAtoms()
    nums = [a.GetAtomicNum() for a in m.GetAtoms() if a.GetAtomicNum() > 1]
    halo = sum(1 for z in nums if z in (9, 17, 35, 53))
    nO = sum(1 for z in nums if z == 8)
    nN = sum(1 for z in nums if z == 7)
    nC = sum(1 for z in nums if z == 6)
    return {
        "MW": Descriptors.MolWt(m),
        "logP": Crippen.MolLogP(m),
        "HBD": Lipinski.NumHDonors(m),
        "HBA": Lipinski.NumHAcceptors(m),
        "RotB": Descriptors.NumRotatableBonds(m),
        "AromRings": rdMolDescriptors.CalcNumAromaticRings(m),
        "Rings": rdMolDescriptors.CalcNumRings(m),
        "TPSA": Descriptors.TPSA(m),
        "FracCSP3": rdMolDescriptors.CalcFractionCSP3(m),
        "heavy": heavy,
        "fracHalo": halo / heavy if heavy else 0,
        "fracO": nO / heavy if heavy else 0,
        "fracC": nC / heavy if heavy else 0,
        "nN": nN,
    }

KEYS = ["MW","logP","HBD","HBA","RotB","AromRings","Rings","TPSA","FracCSP3","heavy","fracHalo","fracO","fracC","nN"]

# reference distributions
refvals = {k: [] for k in KEYS}
for tname, mols in actives.items():
    for mid, smi in mols.items():
        d = desc(smi)
        if d:
            for k in KEYS:
                refvals[k].append(d[k])

print("=== REFERENCE SET descriptor percentiles (n=%d) ===" % len(refvals["MW"]))
pct = {}
for k in KEYS:
    a = np.array(refvals[k])
    p = np.percentile(a, [0.5, 1, 5, 50, 95, 99, 99.5])
    pct[k] = p
    print(f"  {k:10s} p0.5={p[0]:8.2f} p1={p[1]:8.2f} p5={p[2]:8.2f} med={p[3]:8.2f} p95={p[4]:8.2f} p99={p[5]:8.2f} p99.5={p[6]:8.2f}")

id_mols = {**QUERIES, **BACKGROUND}
print("\n=== IN-DOMAIN drugs ===")
print("name".ljust(20)+"".join(k.rjust(9) for k in ["MW","logP","HBD","HBA","RotB","AromR","fracHalo","fracO","nN"]))
id_desc = {}
for n, s in id_mols.items():
    d = desc(s); id_desc[n] = d
    print(n.ljust(20)+f"{d['MW']:9.0f}{d['logP']:9.1f}{d['HBD']:9.0f}{d['HBA']:9.0f}{d['RotB']:9.0f}{d['AromRings']:9.0f}{d['fracHalo']:9.2f}{d['fracO']:9.2f}{d['nN']:9.0f}")

print("\n=== OOD molecules ===")
print("name".ljust(20)+"".join(k.rjust(9) for k in ["MW","logP","HBD","HBA","RotB","AromR","fracHalo","fracO","nN"]))
ood_desc = {}
for n, s in OOD.items():
    d = desc(s); ood_desc[n] = d
    if d is None:
        print(n.ljust(20)+"  unparseable"); continue
    print(n.ljust(20)+f"{d['MW']:9.0f}{d['logP']:9.1f}{d['HBD']:9.0f}{d['HBA']:9.0f}{d['RotB']:9.0f}{d['AromRings']:9.0f}{d['fracHalo']:9.2f}{d['fracO']:9.2f}{d['nN']:9.0f}")

# save percentiles for the rule
json.dump({k: pct[k].tolist() for k in KEYS}, open(os.path.join(HERE,"ref_pctiles.json"),"w"), indent=1)
print("\nWrote ref_pctiles.json")
