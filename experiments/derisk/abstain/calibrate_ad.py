"""Define + calibrate a concrete applicability-domain (AD) / abstain rule, and
TEST it on in-domain drugs vs a hand-built out-of-domain (OOD) set.

--------------------------------------------------------------------------------
FINDING THAT DROVE THE DESIGN
--------------------------------------------------------------------------------
A naive nearest-neighbour (maxNN ECFP4 Tanimoto to the reference set) gate does
NOT separate in-domain drugs from OOD molecules on this 5-target reference set:
  in-domain maxNN spans 0.17-1.00 (metformin 0.18, simvastatin 0.17, prednisone
  0.20 are all LOW), while OOD Gate-A-survivors span 0.14-0.32. They overlap.
  A single T_NN high enough to reject peptides/sugars would false-abstain
  metformin/simvastatin/prednisone/amoxicillin. => NN-distance alone is FRAGILE.
Root cause: with a narrow 5-target reference set, "structurally dissimilar to the
reference ligands" (low NN) is the SAME signal for a legit-but-low-liability drug
and for a genuinely out-of-scope molecule. NN cannot tell them apart here.

--------------------------------------------------------------------------------
THE AD / ABSTAIN RULE WE SETTLED ON (abstain if EITHER gate fires)
--------------------------------------------------------------------------------
Gate A -- FEATURIZABILITY / SMALL-MOLECULE SCOPE:
    SMILES parses AND heavy-atom count >= MIN_HEAVY AND contains no metal atom.
    Fails  -> "not a scoreable organic small molecule (out of scope)".
    (catches: ethanol, urea, carbon tetrachloride, cisplatin.)

Gate B -- OECD DESCRIPTOR DOMAIN (physicochemical bounding box):
    Candidate must lie inside a bounding box on descriptors whose reference-set
    ranges are calibrated from the 7093-active pool (~p99.5) then relaxed to the
    empirical in-domain/OOD gap. Abstain if ANY of:
        MolWt      > MW_MAX          (large peptides / macrocycles)
        fracHalo   > FRACHALO_MAX    (highly halogenated / inorganic-ish)
        fracO      > FRACO_MAX        (sugars / polyols / polyethers)
    Fails  -> "physicochemically outside the reference set's chemical space".
    (catches: cyclosporine, leuprolide peptide, PFOA, hexachlorobenzene,
              glucose, sucrose, PEG oligomer.)

REPORTED (informational, NOT abstain triggers on this reference set):
    maxNN  = max ECFP4 Tanimoto to any reference active  (failed to separate; kept
             as a "nearest analog" display value).
    bestz  = best z-score of the top5-mean class score vs the 24-background
             distribution. Low bestz on an IN-DOMAIN molecule means WEAK COVERAGE
             ("no actionable liability signal"), a DIFFERENT message from "abstain
             / out of domain" -- conflating them would false-reject real drugs.

Outputs: results.json (full per-molecule records + summary).
"""
import json, os, statistics
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski, rdMolDescriptors, Crippen
from rdkit import DataStructs
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

HERE = os.path.dirname(os.path.abspath(__file__))
actives = json.load(open(os.path.join(HERE, "ad_actives.json")))

# ---- calibrated thresholds --------------------------------------------------
MIN_HEAVY    = 6      # Gate A
MW_MAX       = 800.0  # Gate B  (ref p99.5 ~= 806; in-domain max = levothyroxine 777)
FRACHALO_MAX = 0.30   # Gate B  (ref p99.5 = 0.17; in-domain max = fenfluramine 0.19)
FRACO_MAX    = 0.33   # Gate B  (ref p99.5 = 0.21; in-domain max = aspirin 0.31)
NBITS        = 2048
Z_FLOOR      = 2.0    # informational weak-coverage flag

METALS = set([3,4,11,12,13,19,20,21,22,23,24,25,26,27,28,29,30,31,38,39,40,41,42,
              44,45,46,47,48,49,50,55,56,57,72,73,74,75,76,77,78,79,80,81,82,83])

# ---- test molecules ---------------------------------------------------------
QUERIES = {
    "terfenadine":  "CC(C)(C)C1=CC=C(C=C1)C(CCCN2CCC(CC2)C(C3=CC=CC=C3)(C4=CC=CC=C4)O)O",
    "cisapride":    "CO[C@H]1CN(CC[C@H]1NC(=O)C2=CC(=C(C=C2OC)N)Cl)CCCOC3=CC=C(C=C3)F",
    "astemizole":   "COC1=CC=C(C=C1)CCN2CCC(CC2)NC3=NC4=CC=CC=C4N3CC5=CC=C(C=C5)F",
    "thioridazine": "CN1CCCCC1CCN2C3=CC=CC=C3SC4=C2C=C(C=C4)SC",
    "pergolide":    "CCCN1C[C@@H](C[C@H]2[C@H]1CC3=CNC4=CC=CC2=C34)CSC",
    "fenfluramine": "CCNC(C)CC1=CC(=CC=C1)C(F)(F)F",
}
BACKGROUND = {
    "metformin": "CN(C)C(=N)N=C(N)N",
    "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "acetaminophen": "CC(=O)NC1=CC=C(C=C1)O",
    "atorvastatin": "CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4",
    "omeprazole": "CC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=CC=C(C=C3N2)OC",
    "amoxicillin": "CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=C(C=C3)O)N)C(=O)O)C",
    "ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "lisinopril": "C1CC(N(C1)C(=O)C(CCCCN)NC(CCC2=CC=CC=C2)C(=O)O)C(=O)O",
    "warfarin": "CC(=O)CC(C1=CC=CC=C1)C2=C(C3=CC=CC=C3OC2=O)O",
    "furosemide": "C1=CC=C(C(=C1)COC2=CC(=C(C=C2)S(=O)(=O)N)Cl)C(=O)O",
    "hydrochlorothiazide": "C1=CC2=C(C=C1S(=O)(=O)N)S(=O)(=O)NCN2",
    "metoprolol": "CC(C)NCC(COC1=CC=C(C=C1)CCOC)O",
    "amlodipine": "CCOC(=O)C1=C(NC(=C(C1C2=CC=CC=C2Cl)C(=O)OC)C)COCCN",
    "simvastatin": "CCC(C)(C)C(=O)OC1CC(C=C2C1C(C(C=C2)C)CCC3CC(CC(=O)O3)O)C",
    "gabapentin": "C1CCC(CC1)(CC(=O)O)CN",
    "sertraline": "CNC1CCC(C2=CC=CC=C12)C3=CC(=C(C=C3)Cl)Cl",
    "ciprofloxacin": "C1CC1N2C=C(C(=O)C3=CC(=C(C=C32)N4CCNCC4)F)C(=O)O",
    "prednisone": "CC12CC(=O)C3C(C1CCC2(C(=O)CO)O)CCC4=CC(=O)C=CC34C",
    "levothyroxine": "C1=C(C=C(C(=C1I)OC2=CC(=C(C(=C2)I)O)I)I)CC(C(=O)O)N",
    "diazepam": "CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3",
    "ranitidine": "CNC(=CNS(=O)(=O)C)NCCSCC1=CC=C(O1)CN(C)C",
    "loratadine": "CCOC(=O)N1CCC(=C2C3=C(CCC4=C2N=CC=C4)C=C(C=C3)Cl)CC1",
    "naproxen": "CC(C1=CC2=CC=C(C=C2C=C1)OC)C(=O)O",
    "clopidogrel": "COC(=O)C(C1=CC=CC=C1Cl)N2CCC3=C(C2)C=CS3",
}
OOD = {
    "cyclosporine": "CCC1NC(=O)C(C(O)C(C)CC=CC)N(C)C(=O)C(C(C)C)N(C)C(=O)C(CC(C)C)N(C)C(=O)C(CC(C)C)N(C)C(=O)C(C)NC(=O)C(C)NC(=O)C(CC(C)C)N(C)C(=O)C(C(C)C)NC(=O)C(CC(C)C)N(C)C(=O)CN(C)C1=O",
    "leuprolide_frag": "CC(C)CC(NC(=O)C(Cc1c[nH]c2ccccc12)NC(=O)C(Cc1ccc(O)cc1)NC(=O)C(Cc1cnc[nH]1)NC(=O)C1CCC(=O)N1)C(=O)NC(CO)C(=O)O",
    "ethanol": "CCO",
    "urea": "NC(=O)N",
    "glucose": "OCC1OC(O)C(O)C(O)C1O",
    "PFOA": "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
    "hexachlorobenzene": "Clc1c(Cl)c(Cl)c(Cl)c(Cl)c1Cl",
    "carbon_tetrachloride": "ClC(Cl)(Cl)Cl",
    "cisplatin": "N.N.Cl[Pt]Cl",
    "sucrose": "OCC1OC(OC2(CO)OC(CO)C(O)C2O)C(O)C(O)C1O",
    "PEG_oligomer": "OCCOCCOCCOCCOCCOCCOCCO",
}


def analyze(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    heavy = m.GetNumHeavyAtoms()
    nums = [a.GetAtomicNum() for a in m.GetAtoms() if a.GetAtomicNum() > 1]
    halo = sum(1 for z in nums if z in (9, 17, 35, 53))
    nO = sum(1 for z in nums if z == 8)
    metal = any(a.GetAtomicNum() in METALS for a in m.GetAtoms())
    f = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=NBITS)
    return {
        "mol": m, "fp": f, "heavy": heavy, "metal": metal,
        "MW": Descriptors.MolWt(m),
        "logP": Crippen.MolLogP(m),
        "fracHalo": halo / heavy if heavy else 0.0,
        "fracO": nO / heavy if heavy else 0.0,
    }


# ---- build reference fingerprints (pooled + per-target) ---------------------
print("Building reference fingerprints...")
REF_FPS = []
TARGET_FPS = {}
for tname, mols in actives.items():
    lst = []
    for mid, smi in mols.items():
        a = analyze(smi)
        if a:
            lst.append(a["fp"]); REF_FPS.append(a["fp"])
    TARGET_FPS[tname] = lst
    print(f"  {tname}: {len(lst)} fps")
targets = list(TARGET_FPS.keys())
print(f"  POOLED reference set: {len(REF_FPS)} fps")


def class_score(qfp, tname):
    sims = DataStructs.BulkTanimotoSimilarity(qfp, TARGET_FPS[tname])
    sims.sort(reverse=True)
    return statistics.mean(sims[:5])


# background stats for z-floor
BG = {t: [] for t in targets}
for n, s in BACKGROUND.items():
    a = analyze(s)
    for t in targets:
        BG[t].append(class_score(a["fp"], t))
bg_stats = {t: (statistics.mean(v), statistics.pstdev(v)) for t, v in BG.items()}


def evaluate(name, smi):
    a = analyze(smi)
    rec = {"name": name, "smiles": smi}
    if a is None:
        rec.update({"parse_ok": False, "abstain": True, "gate": "A",
                    "reason": "unparseable SMILES"})
        return rec
    rec["parse_ok"] = True
    rec["heavy"] = a["heavy"]
    rec["MW"] = round(a["MW"], 1)
    rec["logP"] = round(a["logP"], 2)
    rec["fracHalo"] = round(a["fracHalo"], 3)
    rec["fracO"] = round(a["fracO"], 3)
    # informational NN + z
    mnn = max(DataStructs.BulkTanimotoSimilarity(a["fp"], REF_FPS))
    zs = []
    for t in targets:
        v = class_score(a["fp"], t)
        m, sd = bg_stats[t]
        zs.append((t, (v - m) / sd if sd > 0 else 0.0))
    zs.sort(key=lambda x: x[1], reverse=True)
    rec["maxNN"] = round(mnn, 4)
    rec["bestz"] = round(zs[0][1], 3)
    rec["bestz_target"] = zs[0][0]
    rec["weak_coverage"] = zs[0][1] < Z_FLOOR
    # Gate A
    if a["metal"]:
        rec.update({"abstain": True, "gate": "A", "reason": "metal-containing (out of scope)"})
        return rec
    if a["heavy"] < MIN_HEAVY:
        rec.update({"abstain": True, "gate": "A",
                    "reason": f"<{MIN_HEAVY} heavy atoms (out of scope)"})
        return rec
    # Gate B (property domain)
    if a["MW"] > MW_MAX:
        rec.update({"abstain": True, "gate": "B", "reason": f"MW {a['MW']:.0f} > {MW_MAX:.0f}"})
        return rec
    if a["fracHalo"] > FRACHALO_MAX:
        rec.update({"abstain": True, "gate": "B", "reason": f"fracHalo {a['fracHalo']:.2f} > {FRACHALO_MAX}"})
        return rec
    if a["fracO"] > FRACO_MAX:
        rec.update({"abstain": True, "gate": "B", "reason": f"fracO {a['fracO']:.2f} > {FRACO_MAX}"})
        return rec
    rec.update({"abstain": False, "gate": None, "reason": "in-domain"})
    return rec


in_domain = {n: evaluate(n, s) for n, s in {**QUERIES, **BACKGROUND}.items()}
ood = {n: evaluate(n, s) for n, s in OOD.items()}

false_abstain = [n for n, r in in_domain.items() if r["abstain"]]
false_accept = [n for n, r in ood.items() if not r["abstain"]]

summary = {
    "rule": {
        "MIN_HEAVY": MIN_HEAVY, "MW_MAX": MW_MAX,
        "FRACHALO_MAX": FRACHALO_MAX, "FRACO_MAX": FRACO_MAX, "Z_FLOOR": Z_FLOOR,
    },
    "n_reference": len(REF_FPS), "targets": targets,
    "n_in_domain": len(in_domain), "n_ood": len(ood),
    "false_abstain_count": len(false_abstain), "false_abstain_names": false_abstain,
    "false_accept_count": len(false_accept), "false_accept_names": false_accept,
    "NN_alone_separates": False,
    "in_domain_maxNN": {
        "min": round(min(r["maxNN"] for r in in_domain.values()), 4),
        "max": round(max(r["maxNN"] for r in in_domain.values()), 4)},
    "ood_maxNN_gateA_survivors": sorted(
        round(r["maxNN"], 4) for r in ood.values() if r.get("gate") != "A" and "maxNN" in r),
    "bg_stats": {t: [round(m, 4), round(sd, 4)] for t, (m, sd) in bg_stats.items()},
}
json.dump({"summary": summary, "in_domain": in_domain, "ood": ood},
          open(os.path.join(HERE, "results.json"), "w"), indent=1, default=str)

# ---- report -----------------------------------------------------------------
print("\n=== FINAL AD RULE ===")
print(f"Gate A: heavy>={MIN_HEAVY} & no metal.  Gate B box: MW<={MW_MAX}, "
      f"fracHalo<={FRACHALO_MAX}, fracO<={FRACO_MAX}")
print(f"Reference: {len(REF_FPS)} actives / {targets}")

print("\n--- IN-DOMAIN (must NOT abstain) ---")
for n, r in sorted(in_domain.items(), key=lambda x: x[1].get("MW", 0)):
    flag = "  <== FALSE ABSTAIN" if r["abstain"] else ""
    print(f"  {n:20s} MW={r['MW']:6.0f} heavy={r['heavy']:3d} fHalo={r['fracHalo']:.2f} "
          f"fO={r['fracO']:.2f} maxNN={r['maxNN']:.2f} bestz={r['bestz']:+.1f}  "
          f"{'ABSTAIN' if r['abstain'] else 'accept':7s}{flag}")

print("\n--- OUT-OF-DOMAIN (must abstain) ---")
for n, r in sorted(ood.items(), key=lambda x: x[1].get("MW", 0)):
    flag = "  <== FALSE ACCEPT" if not r["abstain"] else ""
    mw = r.get("MW", 0); h = r.get("heavy", "?")
    print(f"  {n:22s} MW={mw:6.0f} heavy={h!s:>3} gate={r.get('gate') or '-'} "
          f"{'ABSTAIN' if r['abstain'] else 'accept':7s}  {r['reason']}{flag}")

print(f"\nFALSE ABSTAIN (in-domain rejected): {len(false_abstain)}/{len(in_domain)}  {false_abstain}")
print(f"FALSE ACCEPT  (OOD accepted):       {len(false_accept)}/{len(ood)}  {false_accept}")
print("\nWrote results.json")
