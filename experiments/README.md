# Experiments — reproducible feasibility tests

Backing evidence for the design claims in `../toxicity-assay-recommender.md` §4.
**Environment:** Python 3.13, **RDKit 2026.03.3**, stdlib `urllib` for ChEMBL REST. Conformer seed = 42. Fingerprint = ECFP4 (Morgan r=2, 2048-bit) unless noted.

## What each script does

| Script | Experiment | Key result |
|---|---|---|
| `verify_pairs.py` | Helper: ECFP4 + MACCS Tanimoto + physchem descriptors for a list of pairs | — |
| `run_all.py` | **E1 — 2D similarity** of the 7 verified withdrawn-drug pairs + sanity pairs | Mechanistic pairs ECFP4 **0.05–0.19**; true analogs stay high (pergolide↔cabergoline 0.41, fenfluramine↔norfenfluramine 0.64) → 2D **misses** the mechanistic links |
| `run_pharm.py` | **E2 — feature/pharmacophore 2D** (FCFP4, Gobbi) vs ECFP4 | FCFP lift only **1.3–1.5×**, absolute sim still ~0.15–0.25; Gobbi worse → 2D pharmacophore does **not** rescue |
| `usrcat_experiment.py` | **E3 — 3D shape** (USRCAT, 30 confs ETKDGv3+MMFF), retrieval-rank vs 2D | Partner rank improved **4/7**, worsened **3/7**; one clean rescue (fialuridine–perhexiline 11→1). 3D = weak supplement, **not** a dependable rescue |
| `fetch_actives.py` + `score.py` | **E4 — R4 mechanism linkage** (class-membership scoring vs ChEMBL actives, leave-one-out) | **terfenadine→hERG = 0.57, z = +6.4** vs pairwise 0.19; correct mechanism was the highest-scoring of the 3 tested targets for **4 of 6** drugs, and above the non-binder background for **all 6** → **R4 validated** |
| `exp_fetch_herg.py` + `exp1_specificity.py` | **E5 — specificity** (hERG blockers vs assay-confirmed non-blockers, leave-one-out) | **ROC-AUC 0.894** (blockers median 0.66 vs non-blockers 0.28; 2,494 vs 873); 19/19 common drugs scored low → the score **discriminates**, no false-positive flood |
| `exp2_scaffold.py` | **E6 — scaffold generalization** (leave-one-Bemis-Murcko-cluster-out) | Pooled holdout **AUC 0.913**; congeneric series resilient (0.92–0.98) but a diverse single-ring holdout collapses to **0.667** → real but degrades for truly novel isolated chemotypes |
| `census.py` | **E7a — panel coverage census** (39 Bowes/SAFETYscan targets, # ChEMBL actives pChEMBL≥6) | **37/39 serviceable, 34 rich**; aminergic GPCRs/transporters/enzymes all rich; **BSEP=5, mito=no target, KCNQ1=0** → broad off-target reach, liver is a data desert |
| `phase0_probe.py` + `phase0_potency.py` | **E7b — Phase-0 curation feasibility** (20 tox-failure drugs) | SMILES 20/20 auto; correct mechanism auto-derivable **7/20** (all off-target cases); `failure_reason`/`tier` ~0% structured → hybrid: cheap for off-target drugs, ~2 person-weeks manual for the rest |
| `exp3_multitarget.py` | **E8 — cross-target specificity** (E5 method on non-cardiac targets) | AUC **SERT 0.95, AChE 0.92, MAO-A 0.88** (hERG 0.89 replicates across transporter+enzyme classes); **5-HT2B AUC degenerate** (only 1 ChEMBL inactive) → specificity there via external decoys (12/12 common drugs low). Breadth claim has empirical spine |

`results.json` = E4 output. `exp1_results.json`/`exp2_results.json` = E5/E6. `census_results.csv`/`.json` = E7a (per-target active counts). `phase0_results.json`/`phase0_potency.json` = E7b.

## How to run
```bash
# from a venv with rdkit==2026.03.3
python run_all.py        # E1
python run_pharm.py      # E2
python usrcat_experiment.py   # E3
python fetch_actives.py  # E4 step 1 — pulls ChEMBL actives -> actives.json (network)
python score.py          # E4 step 2 — leave-one-out class scoring -> results.json
python exp_fetch_herg.py # E5/E6 step 1 — hERG blockers + non-blockers from ChEMBL (network)
python exp1_specificity.py   # E5 — specificity AUC
python exp2_scaffold.py      # E6 — scaffold leave-one-cluster-out
python census.py         # E7a — panel coverage census (network)
python phase0_probe.py   # E7b — Phase-0 curation feasibility probe (network)
```
`actives.json` (~217 KB) is regenerable via `fetch_actives.py`, so it's not committed.

## Caveats (see doc §7)
- E4's class-Tanimoto is a **lightweight proxy** for real target-prediction (SEA/QSAR) — inherits ECFP blind spots for small/metabolite-active molecules (fenfluramine miss) and chemotype cross-talk (thioridazine hERG/5-HT2B tie).
- USRCAT is alignment-free global shape; a stricter 3D test = O3A/ROCS overlay (untested).
- ChEMBL fetch depends on pChEMBL≥6 + IC50/Ki/Kd thresholds; cisapride was absent from the hERG pull (fetch incompleteness), which actually strengthens the terfenadine recovery (driven by *other* hERG ligands).
