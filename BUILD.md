# BUILD.md — MVP build spec (for an autonomous Claude loop)

> **You are an autonomous coding agent building the MVP for this project.** This file is your
> single source of truth. Read it fully, then work the milestones in order. It is designed so
> that on every iteration you can: find the first incomplete milestone → implement it → run its
> acceptance check → mark it `[x]` → commit → continue. Stop when M0–M4 all pass and the demo runs.

---

## 0. GOAL (what you are optimizing toward)

**Build a working Streamlit demo that takes a small-molecule SMILES and outputs a *reordered
in-vitro safety-assay plan* — moving the assay most likely to catch a program-ending off-target
liability to the front — with an evidence trail and an LLM-generated med-chemist narrative.**

The product promise, stated as the metric the demo must show on screen:
**"assays-to-culprit" — the killer assay's rank in a default panel vs. in our plan** (e.g. CB1
counter-screen: default rank 15 → ours rank 1).

**Definition of done:** paste any of `rimonabant`, `pergolide`, a novel drug-like molecule, and
`cyclosporine` → the app renders (a) a ranked assay plan with the rank-delta, (b) no-go /
counter-screen / monitor tags, (c) the linked failed drugs + citations, (d) a grounded LLM
narrative, and (e) an abstain message for cyclosporine. Runs locally via `streamlit run`.

---

## 1. CONTEXT — read these first (do not re-derive)

- **The plan & scientific rationale:** [`toxicity-assay-recommender.md`](toxicity-assay-recommender.md) — read §1 (idea), §4 (architecture: R4 engine, Layer-3b severity table, Layer-4 assay map), §6 (demo), §7 (risks).
- **One-pager:** [`SUMMARY.md`](SUMMARY.md).
- **Validated methods & data you REUSE (do not reinvent):**
  - `experiments/score.py` — the R4 scoring backbone: per target, mean-top5 ECFP4 (Morgan r=2, 2048-bit) Tanimoto to that target's ChEMBL actives, z-scored vs a 24-drug BACKGROUND set (the background dict is in this file — copy it).
  - `experiments/derisk/expanded/exp_fetch.py` — the 18-target ChEMBL fetch pattern (pChEMBL≥6, IC50/Ki/Kd).
  - `experiments/derisk/expanded/exp_score.py` — the full 18-target scoring + ranking (adapt this; **remove the leave-one-out** — a real novel candidate removes nothing).
  - `experiments/derisk/drugs.json` — 20 curated failed drugs (smiles, culprit_target, organ, category) → seed for the reference DB.
  - `experiments/derisk/abstain/calibrate_ad.py` — the descriptor-box abstain rule (MW>800 OR halogen-fraction>0.30 OR oxygen-fraction>0.33 OR <6 heavy atoms OR metal → abstain). Reuse it; do NOT use a nearest-neighbour rule (tested, fails).

## 2. HONESTY GUARDRAILS (violating these breaks the pitch — enforce in code + copy)

1. **Never present a score as probability of harm.** Use N×, 0–100 priority, z-score, or rank. No `0.71`-style decimals framed as P(toxic).
2. **The value is buried off-target liabilities, NOT hERG.** hERG is a *validation anchor* only. Lead the demo with **rimonabant (CB1)** and **pergolide (5-HT2B)**. In the UI, when the top hit is hERG, note "standard panels already front-load this."
3. **Frame scope honestly:** the tool addresses *off-target-mediated* failures. It is **blind to metabolite-driven and liver (hepatotox) toxicity** — surface those as `weak-coverage`/`abstain`, never a false negative dressed as "clear."
4. **Show the `known-analog` flag:** if the candidate is highly similar to a known failed drug, say so ("cheap similarity already catches this; our engine adds most for novel chemotypes").
5. **Every reference drug shown must carry a provenance/citation** (from `toxicity-assay-recommender.md` §9 / reference_failures.json).
6. **Metabolite handling:** if a metabolite SMILES is provided, score parent AND metabolite, aggregate by **MAX**.
7. The LLM agent must be **grounded** — pass it only the retrieved structured evidence; instruct it to invent no drugs, mechanisms, or citations.

## 3. STACK & ENVIRONMENT

- **Python deps via `uv`** (NOT pip/conda). Venv at `.venv` (rdkit 2026.03.3 already there). Add packages: `uv pip install --python .venv streamlit anthropic`. Run code with `.venv/bin/python`; run the app with `.venv/bin/streamlit run app/streamlit_app.py`.
- **LLM:** Anthropic Python SDK (`anthropic`), model `claude-opus-4-8` (or `claude-sonnet-5` for speed). Read the `claude-api` skill for current SDK usage. API key from env `ANTHROPIC_API_KEY` — if unset, the agent module must degrade gracefully to a templated report and the app must still run.
- **ChEMBL** REST for panel actives: `https://www.ebi.ac.uk/chembl/api/data/activity.json` (works; see exp_fetch.py).
- **Commit incrementally** (one commit per milestone) on `main`. End commit messages with the Co-Authored-By trailer used elsewhere in this repo.

## 4. TARGET STRUCTURE

```
app/
  data/
    panel.json             # 18 targets: key, chembl_id, assay_name, organ, category, severity, default_action  (author from §4 tables)
    baseline_order.json    # the candidate-independent default panel order (Baseline B, from §4)
    reference_failures.json# curated failed drugs: name, smiles, culprit_target, organ, tier, citation  (seed from drugs.json + §9)
    background.json         # 24-drug background set (copy from experiments/score.py)
    panel_actives.json      # {target_key: {chembl_id: smiles}} — GENERATED by fetch_panel.py (gitignored)
  fetch_panel.py           # builds panel_actives.json from ChEMBL (reuse exp_fetch.py logic)
  build_reference.py       # builds/refreshes reference_failures.json from drugs.json + citations
  core.py                  # score_candidate(smiles, metabolite_smiles=None) -> dict (THE engine)
  agent.py                 # narrative_report(result) -> str  (Anthropic SDK; graceful fallback)
  streamlit_app.py         # the demo UI
  tests/test_core.py       # acceptance checks (runnable with .venv/bin/python -m pytest)
  README.md                # how to run
```
Add `app/data/panel_actives.json` to `.gitignore` (regenerable cache).

## 5. MILESTONES (work in order; each is independently demoable)

### [x] M0 — Data + environment
- `uv pip install --python .venv streamlit anthropic pytest`.
- Author `app/data/panel.json` (18 targets) and `baseline_order.json` from §4's Layer-3b severity table + Layer-4 assay map + the Baseline-B order in `experiments/derisk/expanded/FINDINGS.md`. Panel keys/CHEMBL IDs are listed in exp_fetch.py.
- `build_reference.py` → `reference_failures.json` from `experiments/derisk/drugs.json` + citations in §9.
- Copy the 24-drug BACKGROUND from `experiments/score.py` → `background.json`.
- `fetch_panel.py` → `panel_actives.json` (18 targets, ~600–1000 actives each).
- **Acceptance:** all `app/data/*.json` load; `panel_actives.json` has 18 keys each with >50 SMILES.

### [x] M1 — Core scoring engine (`core.py`)
- `score_candidate(smiles, metabolite_smiles=None)`:
  1. Parse+canonicalize; run the **abstain gate** (reuse calibrate_ad.py rule) → return `{status:"abstain", reason}` if triggered.
  2. Per panel target: mean-top5 ECFP4 Tanimoto to `panel_actives`, z vs `background`. **No LOO.** If metabolite given, MAX-aggregate.
  3. `known-analog` flag: if max Tanimoto to any `reference_failures` drug ≥ 0.5, set it.
  4. Return z-vector + flags.
- **Acceptance (`tests/test_core.py`):** rimonabant → CB1 is top-scoring (z high); pergolide → 5-HT2B in top-3; cyclosporine → `status=="abstain"`; terfenadine → known-analog or hERG top (note as validation anchor).

### [x] M2 — Ranking + reordered plan + evidence (`core.py` cont.)
- Apply Layer-3b severity re-weight → priority per target; map targets→assays (Layer-4); rank → reordered plan.
- Compute **assays-to-culprit delta**: rank of each flagged assay in `baseline_order` vs. ours.
- Attach evidence trail: for each top flagged target, the linked `reference_failures` drugs (name, organ, tier, citation).
- Tag each recommended assay `no-go / counter-screen / monitor` (from severity + known-failure grounding).
- **Acceptance:** `score_candidate("rimonabant SMILES")` returns CB1 counter-screen at rank 1 with default-rank ~15 and an evidence list; output tags present; no P(harm) decimals anywhere.

### [ ] M3 — LLM med-chemist narrative (`agent.py`)
- `narrative_report(result)` → grounded paragraph: names the shared mechanism, the linked failed drugs, and the recommended first assay + why (severity + "resembles drugs withdrawn for X"). Pass ONLY the structured result; instruct no fabrication. Graceful templated fallback if no API key.
- **Acceptance:** for rimonabant, the narrative names CB1 / psychiatric liability, ≥1 real linked failed drug, and recommends the CB1 counter-screen first — with zero invented facts.

### [ ] M4 — Streamlit demo (`streamlit_app.py`)
- SMILES input + example buttons: rimonabant, pergolide, terfenadine (anchor), cyclosporine (abstain), one novel drug-like molecule.
- Render: the reordered assay plan table (assay, our rank, default rank, delta, action tag), the assays-to-culprit headline, the evidence trail, the LLM narrative, and abstain/known-analog banners.
- **Acceptance:** `.venv/bin/streamlit run app/streamlit_app.py` serves; all 5 example molecules render correctly; cyclosporine shows the abstain banner.

### [ ] M5 — Stretch (only if M0–M4 done + time)
- Mechanism-edge network viz (candidate ↔ failed drugs, edges by mechanism, colored by organ).
- Result caching; more reference drugs (scale curation); the assays-to-culprit cost/time weighting (weight assays by real $/turnaround).

## 6. HOW TO OPERATE THIS LOOP (each iteration)
1. Read this file; find the first `[ ]` milestone.
2. Implement it fully (reuse `experiments/` code — do not rewrite the scoring math from scratch).
3. Run its acceptance check with `.venv/bin/python` / pytest / streamlit.
4. If it passes: change `[ ]`→`[x]` here, commit (`M<n>: <what>`), continue to the next.
5. If blocked (e.g. ChEMBL/API down): note it in the milestone, implement the graceful fallback, and continue with mock/cached data rather than stalling.
6. **Stop** when M0–M4 are `[x]` and the demo runs end-to-end. Do not gold-plate; the honest, working vertical slice is the goal.

## 7. DO-NOT list
- Do not rewrite the R4 scoring math — adapt `exp_score.py`.
- Do not add leave-one-out to the live candidate path (LOO was only for retrospective validation).
- Do not claim coverage of hepatotox / metabolite tox — flag as weak-coverage/abstain.
- Do not present scores as probabilities. Do not invent drugs, mechanisms, or citations.
- Do not commit `panel_actives.json` or other regenerable caches (gitignore them).
