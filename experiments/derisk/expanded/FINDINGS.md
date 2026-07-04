# Assay-recovery at n=20 — does the pergolide win GENERALIZE?

**Question.** The n=6 experiment found exactly ONE genuine non-obvious win (pergolide → 5-HT2B, default rank 11 → ours rank 1). Is that a lucky singleton, or does OUR tool systematically surface *buried* (non-cardiac-channel) off-target liabilities that a defensible standard panel misses? Tested on 20 withdrawn/failed drugs: 10 **buried** (culprit = a buried off-target — the hypothesis test) vs 10 **cardiac_herg** (hERG/QT contrast).

## Method (identical scoring backbone to score.py / the n=6 run)
- **Panel:** 18 ChEMBL safety targets → named assays. Actives at pChEMBL≥6, standard_type ∈ {IC50,Ki,Kd}. Fetched counts: hERG 1483, Nav1.5 258, Cav1.2 70, 5-HT2A 1505, 5-HT2B 1214, 5-HT3 740, D2 1523, M2 1397, alpha1A 1354, H1 1285, mu-opioid 1394, SERT 1563, NET 1521, DAT 1387, CB1 1266, AChE 1310, MAO-A 831, COX-1 393. (`exp_actives.json`)
- **Scoring:** per target, mean-top5 ECFP4 (Morgan r=2, 2048-bit) Tanimoto to that target's actives, z-scored vs the SAME 24-drug background as score.py.
- **STRICT leave-one-out:** before scoring a query, its InChIKey **plus every other drug in drugs.json sharing its culprit_target** (mechanistic partners) is removed from EVERY target class. So all recovery comes from *other* ChEMBL ligands, never the reference drugs. (E.g. scoring pergolide removes all 5 5-HT2B drugs from every class.)
- **Three assay rankings, killer-assay rank recorded under each:**
  - **Ours** — 18 assays by candidate z-score, descending.
  - **Baseline B (defensible default, candidate-independent):** hERG, Nav1.5, Cav1.2, 5-HT2A, D2, M2, alpha1A, H1, mu-opioid, SERT, 5-HT2B, 5-HT3, NET, DAT, CB1, AChE, MAO-A, COX-1. Cardiac channels front-loaded (as real panels do); buried off-targets sit low.
  - **Baseline C (base-rate, leave-one-out):** rank assays by culprit frequency across the OTHER 19 drugs (tie-break = B order). **Caveat:** the reference set is enriched for 5-HT2B/hERG by construction, so C structurally over-ranks those two.

## Per-drug result — killer-assay rank (of 18)

| Drug | Category | Culprit | Ours | B | C | killer z |
|---|---|---|---:|---:|---:|---:|
| pergolide | buried | 5-HT2B | **1** | 11 | 2 | +2.25 |
| cabergoline | buried | 5-HT2B | **1** | 11 | 2 | +5.92 |
| fenfluramine | buried | 5-HT2B | 13 | 11 | 2 | +0.35 |
| dexfenfluramine | buried | 5-HT2B | 13 | 11 | 2 | +0.35 |
| methysergide | buried | 5-HT2B | **1** | 11 | 2 | +6.60 |
| rimonabant | buried | CB1 | **1** | 15 | 6 | +10.78 |
| taranabant | buried | CB1 | **1** | 15 | 6 | +9.59 |
| alosetron | buried | 5-HT3 | **1** | 12 | 14 | +2.95 |
| sibutramine | buried | NET | 5 | 13 | 14 | +0.40 |
| mibefradil | buried | Cav1.2 | **1** | 3 | 7 | +10.01 |
| terfenadine | cardiac_herg | hERG | 1 | 1 | 1 | +6.38 |
| cisapride | cardiac_herg | hERG | 2 | 1 | 1 | +5.22 |
| astemizole | cardiac_herg | hERG | 1 | 1 | 1 | +6.22 |
| grepafloxacin | cardiac_herg | hERG | 5 | 1 | 1 | −0.20 |
| sertindole | cardiac_herg | hERG | 1 | 1 | 1 | +9.92 |
| thioridazine | cardiac_herg | hERG | 7 | 1 | 1 | +0.91 |
| terodiline | cardiac_herg | hERG | 1 | 1 | 1 | +1.42 |
| droperidol | cardiac_herg | hERG | 5 | 1 | 1 | +2.91 |
| sparfloxacin | cardiac_herg | hERG | 6 | 1 | 1 | −0.53 |
| mesoridazine | cardiac_herg | hERG | 7 | 1 | 1 | +1.07 |

## Recovery + expected assays-to-culprit

| | Top-1 (Ours/B/C) | Top-3 (Ours/B/C) | Mean assays-to-culprit (Ours/B/C) |
|---|---|---|---|
| **Buried (n=10)** | **7** / 0 / 0 | **7** / 1 / 5 | **3.80** / 11.30 / 5.70 |
| **Cardiac_hERG (n=10)** | 4 / 10 / 10 | 5 / 10 / 10 | 3.60 / 1.00 / 1.00 |

### ★ THE MONEY NUMBER — genuine non-obvious BURIED wins (Ours top-3 AND Baseline B NOT top-3)
**Buried: 6/10** → pergolide, cabergoline, methysergide (5-HT2B valvulopathy); rimonabant, taranabant (CB1 psychiatric); alosetron (5-HT3 colitis).
**Cardiac_hERG: 0/10** (hERG is already rank 1 by default — nothing to add).

The n=6 pergolide singleton is now **6 wins across 3 mechanistically distinct buried liabilities.** These are exactly the assays a cardiac-front-loaded panel buries at ranks 11–15 and that our tool lifts to rank 1.

## Honesty checks — every buried miss, explained (killer NOT in top-3)
- **fenfluramine (rank 13, z=+0.35)** and **dexfenfluramine (rank 13, z=+0.35):** the 5-HT2B agonist is the *metabolite* norfenfluramine; the parent carries no 5-HT2B signal. Known, disclosed limitation (metabolite-active) — R4 cannot see it. Its top hits are hERG/alpha1A/5-HT2A (amphetamine-like), not 5-HT2B.
- **sibutramine (rank 5, z=+0.40):** NET inhibition is driven by its desmethyl metabolites; the parent's strongest signal is DAT (z=+1.55) > NET (z=+0.40). Same metabolite issue + monoamine-transporter cross-talk. Still beats B (rank 13) and C (14), just not into top-3.
- **mibefradil is a hit (Ours rank 1) but NOT a "money" win** because its culprit Cav1.2 is a cardiac ion channel that Baseline B already front-loads to rank 3 — so B also catches it. Correctly excluded from the money count.

**Cardiac_hERG — Ours is WORSE than B (5/10 vs 10/10 top-3), as expected and honest.** For phenothiazines (thioridazine rank 7, mesoridazine 7) the drug's own D2 on-target (z≈+4.8) and other aminergic hits outrank hERG. For fluoroquinolones (grepafloxacin 5, sparfloxacin 6) hERG z is ≈0/negative — structurally distinct antibiotics whose real hERG block isn't captured by ligand similarity. A real panel front-loads hERG at rank 1 regardless, so on cardiac we add nothing and sometimes hurt. The tool's value is NOT hERG.

## Verdict
The crown-jewel claim **generalizes** — but only where the honest framing always said it would: **buried, non-cardiac-channel off-targets.** There, our tool moves the killer assay from a mean rank of **11.3 (default panel) → 3.8**, recovers the culprit in the **top-3 for 7/10** drugs vs **1/10** for the default, and delivers **6 genuine non-obvious wins** spanning three distinct mechanisms (5-HT2B valvulopathy, CB1 psychiatric, 5-HT3 colitis). Baseline C, even while structurally over-ranking 5-HT2B by construction, manages only 5/10 and is helpless on CB1/5-HT3/NET/Cav1.2. On cardiac/hERG the tool adds nothing (hERG is default rank 1) and is sometimes worse — stated plainly, not hidden. The two clean misses (fenfluramine/dexfenfluramine) and one partial (sibutramine) are all metabolite-active drugs, the documented boundary of a parent-structure method.

## Files
- `exp_fetch.py` — ChEMBL fetch for 18 targets (→ `exp_actives.json`)
- `exp_score.py` — strict-LOO scoring + three-way ranking, per-category recovery, money number (→ `results.json`)
- `results.json` — full z-vectors, per-drug ranks (Ours/B/C), recovery counts, mean expected-assays-to-culprit, money drugs
