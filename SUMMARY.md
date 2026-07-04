# Toxicity → Assay Recommender — one-page summary

*Full detail: [`toxicity-assay-recommender.md`](toxicity-assay-recommender.md) · evidence: [`experiments/`](experiments/)*

## What it is
Paste a drug candidate's structure → the tool **re-orders your in-vitro safety-assay plan so the test most likely to catch a real liability runs first.** It does this by linking the candidate to **drugs that actually failed in the clinic for toxicity**, via the off-target mechanism they share. It prioritizes experiments — it does **not** claim a molecule is safe.

## How it works (3 layers, plain version)
1. **Find the mechanism, not the lookalike.** A candidate rarely *looks* like the drug it will fail alongside (2D similarity of real failure-pairs is near-zero, 0.05–0.19). Instead we score the candidate against the **whole class of known ligands** for each "safety panel" off-target (hERG, 5-HT2B, …). Resembling *some* hERG blocker ⇒ likely hERG liability — even when it resembles no single one.
2. **Weight by consequence, not just likelihood.** A hit isn't automatically a killer. Each off-target carries a **severity/action tag** (no-go / counter-screen / monitor), boosted when the candidate resembles drugs that were *actually withdrawn* for that mechanism.
3. **Recommend the assay + tell the truth about confidence.** Map mechanism → named assay, tagged **mechanism-linked / model-predicted / abstain**. A validated engine (M1) covers off-target liabilities across cardiac/CNS/GI/endocrine/immune; wrapped models (M2) cover the liver/mito gap at a lower, clearly-labeled tier.

## What we already proved (cheap experiments, all in `experiments/`)
- **It discriminates:** separates toxic from non-toxic — **ROC-AUC 0.89 (hERG)**, replicated **SERT 0.95, AChE 0.92, MAO-A 0.88** across transporter/enzyme classes. No "flags everything" failure.
- **It's broad:** **37/39** safety-panel targets have enough data to serve.
- **It finds the non-obvious:** recovers 7 real withdrawn-drug pairs that 2D structure search completely misses (e.g. terfenadine↔cisapride, z = +6.4 to the hERG class, leave-one-out).
- **We know the limits:** liver (BSEP/mito) is a data desert for this method; truly novel scaffolds degrade (AUC → 0.67).

## Main selling points
- **Grounded in real clinical failures**, not textbook alerts — "you resemble drugs that *died* for this, run this assay first."
- **A number, not a vibe:** the demo shows we move the killer assay from rank N → **top-3**.
- **Honest by construction:** evidence tiers + an **abstain** option turn the coverage limit into a trust feature.
- **Broad but validated:** works across many off-target liabilities, with discrimination measured on 4 targets — not a single-endpoint predictor.
- **Cheap to build:** retrieval + light reasoning over public data (ChEMBL, DILIrank/DICTrank); no model training needed for v1.

## Main risks (and the honest answer)
- **Curation is the real bottleneck** — labeling *why* a drug failed (organ/mechanism/provenance) is ~2 person-weeks for a 200-drug set, not a same-day merge. → scope to the auto-derivable off-target drugs first.
- **Liver/mito/reactive-metabolite ≠ off-target binding** — the core engine can't reach them. → served by lower-tier wrapped models (M2), clearly flagged, not oversold.
- **Novel chemotypes degrade** — it's a similarity method; a molecule far from all known ligands gets weak signal. → abstain + confidence flags.
- **Breadth partly unmeasured** — 4 targets validated, ~33 rich ones still assumed. → Phase-3 work; don't overclaim "all 37."
- **Some GPCRs lack clean negatives** in ChEMBL (e.g. 5-HT2B) → specificity there rests on external decoys, not a clean AUC.

## Status
De-risking **done** (verdict: viable, broad within off-target space, honest liver/mito gap). No product code yet. Build order: **M1 off-target core + assay-ranking demo → M2 breadth → hardening.**
