# Toxicity Risk → Assay Recommender (MVP)

Broader orientation: [`../README.md`](../README.md).

Paste a small-molecule SMILES → a **reordered in-vitro safety-assay plan** that moves the
assay most likely to catch a program-ending **off-target** liability to the front, with an
evidence trail (linked failed drugs + citations) and a grounded LLM med-chemist narrative.

Product metric on screen — **assays until culprit**: the killer assay's rank in a default panel
vs. in our plan (e.g. rimonabant CB1 counter-screen: default rank 15 → ours rank 1).

## Run

```bash
# deps (uv, not pip/conda)
uv pip install --python .venv streamlit anthropic pytest      # rdkit already in .venv

# build data (one-time; panel_actives.json is a gitignored ChEMBL cache)
.venv/bin/python app/build_reference.py      # -> app/data/reference_failures.json
.venv/bin/python app/fetch_panel.py          # -> app/data/panel_actives.json (18 targets)

# run the demo
.venv/bin/streamlit run app/streamlit_app.py

# tests
.venv/bin/python -m pytest app/tests -q
```

Set `ANTHROPIC_API_KEY` for the LLM narrative (model `claude-opus-4-8`); without it the app
degrades to a deterministic, fully-grounded templated report and still runs.

## What to try

- **Rimonabant** — CB1 counter-screen jumps to rank 1 (default 15). Buried liability; the win.
- **Pergolide** — 5-HT2B counter-screen surfaced (default 11). Valvulopathy.
- **Terfenadine** — hERG anchor; app notes "standard panels already front-load this".
- **Novel aryl-piperazine** — in-domain, *not* a known analog; the novel-chemotype regime.
- **Cyclosporine** — abstain (MW 1203, outside applicability domain).

## Method (reused from `experiments/`, not reinvented)

The main engine is **off-target class matching**. The app ships an 18-target MVP subset of
the 37 serviceable Bowes/SAFETYscan panel targets (coverage census in `experiments/`).
Per panel target: mean-top5 ECFP4 (Morgan r=2, 2048-bit) Tanimoto to that target's ChEMBL
actives, z-scored vs a 25-drug background (`experiments/score.py`). **No leave-one-out** on the
live path; leave-one-out exists only for demo/validation mode. Layer-3b severity re-weight →
priority; Layer-4 target→assay map; descriptor-box abstain gate
(`experiments/derisk/abstain/calibrate_ad.py`); known analog flag at Tanimoto ≥ 0.5; metabolite
MAX-aggregation. Liver/mito/reactive-metabolite checks are lower-confidence organ-tox modules and
do not change the validated off-target headline.

## Honesty guardrails (enforced in code + copy)

- Scores are **z-score / N× / rank / 0–100 priority — never P(harm)**.
- Value is **buried off-targets**, not hERG (a validation anchor only).
- Scope: **off-target-mediated** failures; **blind to metabolite-driven and hepatotoxicity**
  — surfaced as weak-coverage / abstain, never a false "clear".
- Every reference drug shown carries a **provenance/citation**.
- The LLM agent is **grounded** — it sees only retrieved structured evidence and is instructed
  to invent no drugs, mechanisms, or citations.

## Files

```
app/
  data/            panel.json, baseline_order.json, reference_failures.json, background.json,
                   panel_actives.json (gitignored cache)
  fetch_panel.py   builds panel_actives.json from ChEMBL
  build_reference.py  builds reference_failures.json from experiments/derisk/drugs.json + citations
  build_outcomes.py   validates/normalizes outcome_actives.json + reactive_alerts.json
  core.py          score_candidate() + build_plan() — the engine
  outcome_modules.py lower-confidence liver/mito/reactive-metabolite checks
  render.py        molecule drawings + mechanism graph
  validation.py    the measured-numbers panel (cited to experiments/derisk FINDINGS)
  agent.py         narrative_report() — grounded Anthropic SDK, templated fallback
  streamlit_app.py the demo UI
  tests/          core, outcome-module, and Streamlit render checks
```
