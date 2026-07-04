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
| `fetch_actives.py` + `score.py` | **E4 — R4 mechanism linkage** (class-membership scoring vs ChEMBL actives, leave-one-out) | **terfenadine→hERG = 0.57, z = +6.4** vs pairwise 0.19; true mechanism top-ranked **4/6**, above background **6/6** → **R4 validated** |

`results.json` = E4 output (per-drug class scores, background stats, actives counts: hERG 1,483 / 5-HT2B 1,214 / Nav1.5 258; BSEP only 5 → unusable, fell back).

## How to run
```bash
# from a venv with rdkit==2026.03.3
python run_all.py        # E1
python run_pharm.py      # E2
python usrcat_experiment.py   # E3
python fetch_actives.py  # E4 step 1 — pulls ChEMBL actives -> actives.json (network)
python score.py          # E4 step 2 — leave-one-out class scoring -> results.json
```
`actives.json` (~217 KB) is regenerable via `fetch_actives.py`, so it's not committed.

## Caveats (see doc §7)
- E4's class-Tanimoto is a **lightweight proxy** for real target-prediction (SEA/QSAR) — inherits ECFP blind spots for small/metabolite-active molecules (fenfluramine miss) and chemotype cross-talk (thioridazine hERG/5-HT2B tie).
- USRCAT is alignment-free global shape; a stricter 3D test = O3A/ROCS overlay (untested).
- ChEMBL fetch depends on pChEMBL≥6 + IC50/Ki/Kd thresholds; cisapride was absent from the hERG pull (fetch incompleteness), which actually strengthens the terfenadine recovery (driven by *other* hERG ligands).
