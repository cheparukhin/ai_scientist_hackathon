# ProTox-3.0 — complete predicted-endpoint inventory

Source (3 concordant): live Models-page dropdown (Wayback snapshot 2024-05 / 2026 of tox.charite.de/protox3/index.php?site=models),
the official API script `ALL_MODELS` string, and the NAR-2024 paper (52:W513, PMC11223834). Captured 2026-07-04.

## 1. Acute toxicity (1)
- Oral toxicity: predicted LD50 (mg/kg) + toxicity class 1-6.

## 2. Organ / endpoint models (13, binary Active/Inactive + confidence)
Hepatotoxicity (dili) · Neurotoxicity (neuro) · Nephrotoxicity (nephro) · Respiratory tox (respi) ·
**Cardiotoxicity (cardio)** · Carcinogenicity · Immunotoxicity · Mutagenicity · Cytotoxicity ·
BBB-barrier · Ecotoxicity · Clinical toxicity · Nutritional toxicity.
  - Per NAR-2024, the **cardiotoxicity** model is trained on hERG-blocker data (so hERG liability surfaces here as a coarse organ flag).

## 3. Tox21 pathways (12, nuclear-receptor + stress-response)
NR: AhR, AR, AR-LBD, Aromatase, ER, ER-LBD, PPAR-gamma.
SR: nrf2/ARE, HSE, MMP (mito membrane potential), p53, ATAD5.

## 4. Molecular Initiating Events — MIE (14)
THR-alpha, THR-beta, TTR, RYR (ryanodine), GABAR, NMDAR, AMPAR, KAR (kainate),
AChE, CAR, PXR, NADHOX, **VGSC (voltage-gated Na+ channel / Nav)**, NIS (Na/I symporter).

## 5. Metabolism CYPs (6)
CYP1A2, CYP2C19, CYP2C9, CYP2D6, CYP3A4, CYP2E1.

## 6. Toxicity Targets (15, pharmacophore "toxicophore" fit, class 0-3)
Definition (ProTox-II NAR-2018, PMC6031011): "all targets belonging to the Novartis in-vitro safety panel...
from 73 targets, **15 were found to have at least one experimental structure representing a human
protein-ligand complex**." i.e. a **solved-crystal-structure filter**.
Confirmed members surfaced in the literature: **Adrenergic beta-2 receptor** (7 toxicophores),
**Dopamine D3 receptor** (1 toxicophore). Full 15 not enumerated on any accessible page/paper.
Output = average pharmacophore fit + ligand similarity, reported as a binding-probability CLASS 0-3.

## What is ABSENT as a dedicated target/endpoint
- No standalone **hERG/KCNH2** target output (only inside the cardio ML model; no hERG crystal structure until 2017 cryo-EM, after target selection).
- No **5-HT3** (cys-loop pentamer; only non-human structures existed at selection time).
- No **NET / SERT / DAT** monoamine-transporter output (no human transporter structure at selection time).
- No **Cav** L-type calcium channel target (RYR + VGSC are the only channels, as MIEs).
- **5-HT2B** and **CB1**: human structures exist (5-HT2B 4IB4 2013; CB1 2016), so either COULD be among the
  unenumerated 15 — cannot be positively confirmed or excluded from documentation alone.

## Output format (all sections)
Flat table: input | type | Target | Prediction | Probability. No ranking, no assay recommendation,
no linkage to clinically-failed drugs, no severity/actionability tier, no reordering.
