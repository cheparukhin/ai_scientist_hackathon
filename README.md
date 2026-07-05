# Toxicity Risk Assay Recommender

This repo contains a working Streamlit demo plus the research evidence behind it.

Paste a small-molecule SMILES and the app reorders an in-vitro safety panel so the
test most likely to catch a program-ending off-target liability runs first. It is
a test-prioritization tool, not a safety classifier.

## Quick Start

```bash
# install app dependencies into the existing uv-managed venv
uv pip install --python .venv streamlit anthropic pytest

# build or refresh local data
.venv/bin/python app/build_reference.py
.venv/bin/python app/fetch_panel.py

# run the demo
.venv/bin/streamlit run app/streamlit_app.py

# run the test suite
.venv/bin/python -m pytest app/tests -q
```

`app/data/panel_actives.json` is a generated ChEMBL cache and is intentionally
gitignored. Set `ANTHROPIC_API_KEY` to enable the LLM narrative; without it the
app uses a deterministic grounded fallback.

## Read This First

| File | Use it for |
|---|---|
| [`SUMMARY.md`](SUMMARY.md) | One-page product and evidence summary. |
| [`app/README.md`](app/README.md) | How to run, test, and understand the demo app. |
| [`experiments/README.md`](experiments/README.md) | Map of the reproducible experiments and their outputs. |
| [`toxicity-assay-recommender.md`](toxicity-assay-recommender.md) | Long-form research memo and pitch rationale. |
| [`BUILD.md`](BUILD.md), [`BUILD2.md`](BUILD2.md) | Historical implementation specs for the MVP and phase-2 polish. |

## How The Main Engine Works

The main engine is **off-target class matching**:

1. For each safety-panel target, collect known active molecules from ChEMBL.
2. Compare the candidate to the whole active class, not just one failed drug.
3. Use the mean top-5 ECFP4 Tanimoto similarity as the target score.
4. Convert that score to a z-score against a fixed 25-drug background of ordinary marketed drugs.
5. Rank assays by that z weighted by a per-target severity tier, then map the top target to its named assay.

### Worked example — rimonabant at CB1 (the exact numbers the app shows)

Rimonabant was withdrawn for psychiatric/suicidality acting through the **CB1** receptor —
a liability a default panel wouldn't counter-screen until position **#15**. Every step below
is reproducible under *Show the calculation* in the app:

| Step | Value |
|---|---|
| Candidate fingerprint | ECFP4 (Morgan r=2, 2048-bit) |
| Known CB1 binders compared (ChEMBL, pChEMBL ≥ 6) | 1,266 |
| Mean of the top-5 Tanimoto to that class (raw score) | 0.940 |
| Background over 25 ordinary drugs: mean / SD | 0.229 / 0.062 |
| **z = (0.940 − 0.229) / 0.062** | **+11.4** |
| Flagged (z ≥ 2)? | yes |
| Default panel rank → our rank | **#15 → #1** |

In words: rimonabant is **~11 standard deviations** more similar to known CB1 binders than a
typical marketed drug. That is **similarity-enrichment, not a probability of harm** — no dose,
exposure, or metabolism is modelled. Turn on **Demo mode** and the tool hides rimonabant *and its
withdrawn cousins*, yet still recovers CB1 at **z ≈ +10.8** from *other* CB1 ligands — evidence the
signal is the shared mechanism, not memorization of the drug itself.

The early research compared this against 2D fingerprint matching,
feature/pharmacophore fingerprints, and 3D shape matching. Off-target class
matching is the useful core for novel chemotypes because it can recover a shared
binding mechanism even when no single failed drug looks similar.

## What The Evidence Supports

- For buried off-target liabilities, the retrospective validation moved the
  culprit assay from mean rank **11.3** in a default panel to **3.8** in our
  plan.
- The strongest use case is a **novel chemotype** whose risk is not obvious from
  nearest-neighbor similarity to a known failed drug.
- The app intentionally says when it has low confidence: `abstain`,
  known analog, weak coverage, and lower-confidence organ-tox modules are
  part of the product, not afterthoughts.

## Important Bounds

- The main validated claim is for **off-target-mediated failures**, not every
  toxicity failure.
- hERG is a validation anchor, but not a product win: standard panels already
  run hERG early.
- Liver, mitochondrial, and reactive-metabolite findings are lower-confidence
  organ-tox modules. They are hypotheses to confirm experimentally, not
  probabilities of harm.
