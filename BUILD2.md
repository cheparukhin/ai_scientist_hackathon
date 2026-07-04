# BUILD2.md — Phase-2 build spec (demo polish + lower-confidence outcome breadth)

> **You are an autonomous coding agent extending a working MVP.** The MVP (BUILD.md M0–M4)
> is complete: `app/core.py` (`score_candidate`, `build_plan`), `app/agent.py`
> (`narrative_report`), `app/streamlit_app.py`, `app/tests/test_core.py` (10 passing).
> This file is your single source of truth for the NEXT phase. Read it fully, then work
> milestones **N1→N7 in order**. On each iteration: find the first incomplete `[ ]` →
> implement it → run its acceptance check → mark it `[x]` → commit → continue. Stop when
> N1–N6 are `[x]` (N7 is stretch) and `.venv/bin/python -m pytest app/tests -q` is green.

---

## 0. GOAL

Close the gap between the working MVP and the **§6 investor demo** in
`toxicity-assay-recommender.md`, and add lower-confidence outcome breadth for the
metabolism/liver/mito gap the validated off-target core is blind to.

Two workstreams:
- **Demo polish (N1–N4, N7):** structure images, a leave-one-out "demo mode" that turns the
  rimonabant/pergolide screens into the intended *value moment*, an in-app validation panel
  with the measured numbers, a scope-&-limits expander, and a mechanism-edge network viz.
- **Outcome breadth (N5–N6):** reactive-metabolite structural alerts + hepatotox/mito read-across
  modules, fused into the plan as `model-predicted` (tier-2) rows — never merged into the
  validated off-target headline.

**Definition of done:** N1–N6 `[x]`; `app/tests` all green; `.venv/bin/streamlit run
app/streamlit_app.py` renders, for rimonabant: candidate structure image, the assays-until-culprit
headline, a **leave-one-out toggle** that (when on) recovers CB1 with the known analog flag OFF,
a validation panel, a scope expander, and a distinct **model-predicted organ-tox** section;
and a hepatotoxicant (troglitazone) surfaces a hepatotox + mito model-predicted flag while its
validated off-target headline is unchanged.

---

## 1. CONTEXT — read these first (do not re-derive)

- **Plan:** `toxicity-assay-recommender.md` — §5 Phase-2 (outcome breadth), §6 demo (the
  leave-one-out framing + trust signals), §7 risks (esp. lower-confidence organ-tox modules and
  narrow population), §4 severity/actionability and mechanism-to-assay mapping.
- **Validated numbers to display (do NOT invent — cite these files):**
  `experiments/derisk/expanded/FINDINGS.md` — buried mean assays until culprit **11.30 → 3.80**,
  top-3 **7/10 vs 1/10**, **6/10** money wins (pergolide, cabergoline, methysergide → 5-HT2B;
  rimonabant, taranabant → CB1; alosetron → 5-HT3). Per-drug Ours/B/C ranks are in that file's
  table. Per-target AUC / scaffold-split / ablation are in `toxicity-assay-recommender.md` §4
  (hERG 0.89, SERT 0.95, AChE 0.92, MAO-A 0.88; scaffold 0.913 pooled / 0.667 novel-isolated;
  ablation naive 0/20 vs off-target class matching 12/20).
- **Leave-one-out logic to reuse (do NOT rewrite):** `experiments/derisk/expanded/exp_score.py`
  (`loo_exclude`, `by_culprit`) — exclude the query's InChIKey **plus every reference drug
  sharing its `culprit_target`** from every target class. This is the *demo-only* path.
- **Existing engine internals:** `app/core.py` `_Engine` stores `target_fps[t] = [fp, ...]`
  (no InChIKeys yet — N2 adds them), `bg_stats`, `ref` (name, fp, meta). `build_plan(result)`
  returns `{rows, headline, flagged}`.

## 2. HONESTY GUARDRAILS (carry over from BUILD.md §2 + outcome-module-specific — enforce in code + copy)

1. **Never present a score as probability of harm** — z / N× / 0–100 priority / rank only.
   This applies to liver/mito read-across too (it is similarity-enrichment, not P(harm)).
2. **Outcome modules are a lower evidence tier than the off-target core.** Tag every organ-tox item
   `evidence_tier="model-predicted"` and render it visually separated and below the mechanism-linked
   off-target plan. Surface the caveat "less validated than the off-target core" wherever these
   modules appear. **Never merge them into the assays-until-culprit headline** — that metric belongs
   only to the validated off-target claim.
3. **Reactive-metabolite output is an *alert*, not a prediction** — a structural liability
   hypothesis that needs experimental confirmation. Say so.
4. **Every curated outcome drug carries a provenance/citation** (as reference_failures does).
5. **Leave-one-out demo mode must not change the live/novel path.** A real novel candidate removes
   nothing; leave-one-out only fires when the user explicitly enables demo mode AND the input matches a
   known reference drug. Default is OFF.
6. **Keep the buried > hERG framing and the known analog honesty** (BUILD.md §2.2–2.4).

## 3. STACK & ENVIRONMENT

- Deps via `uv` (`.venv`). Everything needed is already installed: `rdkit` (incl.
  `rdkit.Chem.FilterCatalog` and `Draw`), `streamlit` (`st.image`, `st.graphviz_chart`,
  `st.expander`, `st.bar_chart`), `pandas`, `anthropic`, `pytest`. **Do not add new
  dependencies** unless a milestone explicitly authorizes it (none do — TDC/ProTox are
  optional stretches only, behind graceful guards).
- Run code with `.venv/bin/python`; app with `.venv/bin/streamlit run app/streamlit_app.py`;
  tests with `.venv/bin/python -m pytest app/tests -q`.
- Verify UI milestones headlessly with `streamlit.testing.v1.AppTest` (no browser needed).
- Commit one commit per milestone on `main`. End messages with the Co-Authored-By trailer.

## 4. TARGET STRUCTURE (new/changed files)

```
app/
  render.py               # NEW  mol_png(smiles) + mechanism_graph_dot(result, plan)  (viz helpers)
  outcome_modules.py      # NEW  lower-confidence reactive-metabolite alerts + hepatotox/mito read-across
  data/
    outcome_actives.json  # NEW  curated DILI+ / mito-toxicant sets (name, smiles, endpoint, citation) — COMMITTED
    reactive_alerts.json  # NEW  curated bioactivation SMARTS (name, smarts, note, citation)      — COMMITTED
  build_outcomes.py       # NEW  validates/normalizes outcome_actives.json + reactive_alerts.json (SMILES/SMARTS parse)
  core.py                 # CHANGED  add leave-one-out support (InChIKey exclusion) to _Engine + score_candidate
  agent.py                # CHANGED  tier-labelled mention of organ-tox findings (grounded)
  streamlit_app.py        # CHANGED  structure images, leave-one-out toggle, validation panel, scope expander,
                          #          organ-tox section, network viz
  tests/
    test_core.py          # CHANGED  keep 10 green; add N2 leave-one-out tests
    test_outcomes.py      # NEW  outcome-module acceptance
    test_app.py           # NEW  AppTest render checks (structure/leave-one-out/validation/outcome sections)
```

---

## 5. MILESTONES (work in order; each independently demoable; commit per milestone)

### [x] N1 — 2D structure depiction (`render.py`)
- `render.py`: `mol_png(smiles, size=(320,240)) -> bytes|None` using RDKit
  `rdMolDraw2D.MolDraw2DCairo` (fallback `Draw.MolToImage` → PNG bytes); return `None` on
  unparseable SMILES.
- In `streamlit_app.py`: show the candidate structure (`st.image`) near the top; for an OK
  result also show the **headline linked failed drug's** structure side-by-side (look up its
  SMILES from `reference_failures.json`), captioned with the drug name. Never crash on a bad
  SMILES (guard `None`).
- **Acceptance (`tests/test_app.py`):** AppTest on rimonabant renders ≥1 `st.image` element
  (candidate structure); `mol_png(EXAMPLES["rimonabant"])` returns non-empty bytes;
  `mol_png("not_a_smiles")` returns `None`.

### [x] N2 — Leave-one-out "demo mode" (`core.py`)
- Extend `_Engine`: store `target_fps[t]` as `[(inchikey, fp), ...]` (compute InChIKey via
  `Chem.MolToInchiKey`). Keep a fast path when no exclusion set is given (unchanged behaviour).
- `_class_score(qfp, target, exclude_iks=frozenset())` and `score_fp(fp, exclude_iks=...)` skip
  fps whose InChIKey ∈ `exclude_iks`. **Recompute `bg_stats` is NOT needed** — background is
  candidate-independent; keep it as-is (leave-one-out only excludes from the query's per-target scores and
  from the nearest-ref computation).
- `score_candidate(smiles, metabolite_smiles=None, loo=False)`:
  - `loo=False` → identical to today (live/novel path; removes nothing).
  - `loo=True` → build the exclusion set from the candidate: its own InChIKey **plus** every
    `reference_failures` drug sharing the candidate's matched culprit target (mirror
    `exp_score.loo_exclude`). Match the candidate to a reference drug by InChIKey; if it isn't a
    known reference drug, leave-one-out exclusion is just its own InChIKey. Exclude these IKs from every
    target's `_class_score` AND from `nearest_ref` (so a self-match can't set known analog).
- UI: a checkbox **"Demo mode — leave-one-out (score this known drug as if novel: remove it and
  its mechanistic partners from the DB)"**, default OFF, passed through to `score_candidate`.
  When ON and the input is a known reference drug, show a small badge explaining it.
- **Acceptance (`tests/test_core.py`):** `score_candidate(EXAMPLES["rimonabant"], loo=True)` →
  `best_target=="CB1_CHEMBL218"`, CB1 `z >= 6.0`, `flags["known_analog"] is False`,
  `nearest_analog["name"] != "rimonabant"`. And `loo=False` still gives `known_analog True`
  (regression). Existing 10 tests stay green.

### [x] N3 — Scope & limits expander (`streamlit_app.py`)
- Add `st.expander("Scope & limits — read before trusting this")` containing the §7/§0 caveats
  in plain language: (a) covers *off-target-mediated* failures only; blind to metabolite-driven
  and liver tox → those are the *model-predicted* outcome modules at a lower tier; (b) narrow
  addressable population (most withdrawals are off-panel hepato/on-target — §7.10); (c) no value
  on hERG (validation anchor only); (d) advantage is conditional on novelty (known analog flag).
- **Acceptance:** AppTest → an expander labelled with "Scope" exists (`at.expander` non-empty).

### [x] N4 — Validation panel (`streamlit_app.py`)
- Add a "Validation — measured on 20 historical failures (leave-one-out)" panel/expander with
  the FINDINGS numbers **verbatim, cited to `experiments/derisk/expanded/FINDINGS.md`**:
  buried mean assays until culprit **11.3 → 3.8**, top-3 **7/10 vs 1/10 (default)**, **6 non-obvious
  wins**; per-target AUC (hERG 0.89, SERT 0.95, AChE 0.92, MAO-A 0.88); scaffold-split 0.913
  pooled / 0.667 novel-isolated; ablation naive 0/20 vs off-target class matching 12/20. Include a small `st.bar_chart`
  of default-panel rank vs our rank for the 10 buried drugs (hard-code the rank pairs from the
  FINDINGS table — pergolide 11→1, cabergoline 11→1, fenfluramine 11→13, dexfenfluramine 11→13,
  methysergide 11→1, rimonabant 15→1, taranabant 15→1, alosetron 12→1, sibutramine 13→5,
  mibefradil 3→1). No fabricated numbers.
- **Acceptance:** AppTest → panel present and contains the strings "11.3" and "3.8"; the chart
  data has 10 rows. A test asserts the hard-coded rank pairs match the FINDINGS.md table
  (parse the table, compare).

### [x] N5 — Lower-confidence outcome modules (`outcome_modules.py` + curated data)
- **`data/reactive_alerts.json`** (committed): curated bioactivation SMARTS with citations
  (Kalgutkar 2005; Stepan 2011; Claydon/Park DILI alerts). Cover at least: aromatic nitro,
  aromatic amine/anilide, thiophene, furan, hydrazine/hydrazide, thiourea, epoxide,
  Michael-acceptor enone, quinone / para-aminophenol (quinone-imine former), α-halo carbonyl,
  terminal aryl alkyne, free thiol. Each: `{name, smarts, note, citation}`.
- **`data/outcome_actives.json`** (committed): read-across reference sets, each drug
  `{name, smiles, endpoint, tier, citation}`:
  - `hepatotox` (DILI-positive, ~20 drugs, organ-concern tier, cite DILIrank / FDA-LTKB): troglitazone,
    nefazodone, trovafloxacin, tolcapone, ketoconazole, diclofenac, bromfenac, isoniazid,
    valproic acid, amiodarone, felbamate, ticlopidine, pemoline, nimesulide, flutamide,
    benoxaprofen, methotrexate, dantrolene, labetalol, ebrotidine.
  - `mito` (~12 drugs, cite Nadanaciva/Will 2007; Aleo 2014): troglitazone, nefazodone,
    tolcapone, perhexiline, flutamide, amiodarone, tamoxifen, benzbromarone, cerivastatin,
    valproic acid, chlorpromazine, felbamate.
  - **Verify every SMILES parses** (build_outcomes.py); fix any that don't from a trusted source.
- **`build_outcomes.py`**: load both JSONs, assert all SMILES/SMARTS parse in RDKit, report counts.
- **`outcome_modules.py`**:
  - `reactive_alerts(smiles) -> [{name, note, citation}]` — RDKit `FilterCatalog`
    (BRENK+NIH) matches **plus** the curated SMARTS matches; dedupe by name.
  - `outcome_scores(smiles) -> {endpoint: {z, raw, nearest:{name,sim,citation}, tier,
    flagged}}` for `hepatotox` and `mito` — reuse the off-target class-matching backbone (mean-top5 ECFP4 Tanimoto to
    that endpoint's actives, z vs the SAME 24-drug background as core). `flagged = z >= 2.0`.
  - `outcome_panel(smiles) -> {reactive: [...], endpoints: {...}, tier:"model-predicted"}`.
    All z-based, no probabilities.
- **Acceptance (`tests/test_outcomes.py`):** troglitazone → `hepatotox` flagged AND `mito`
  flagged; a nitroaromatic (e.g. `"O=[N+]([O-])c1ccc(N)cc1"`) → `reactive_alerts` non-empty and
  includes a nitro/aniline alert; metformin → hepatotox NOT flagged and no reactive alert;
  every endpoint dict has `tier=="model-predicted"` and no key named `probability`/`p_toxic`.

### [x] N6 — Fuse outcome modules into UI + narrative
- `streamlit_app.py`: below the validated off-target plan/evidence, add a clearly-separated section
  **"Metabolism / organ-tox — model-predicted (lower evidence tier, less validated than the
  off-target core)"** showing: reactive-metabolite alerts (as chips/list with the matched
  liability + citation + "hypothesis, confirm experimentally"), and hepatotox/mito rows
  (endpoint, z, nearest DILI/mito analog + citation, `flagged` tag). Show even when nothing is
  flagged (state "no model-predicted liability flagged — still not a clean bill of health").
  **Do not change the assays-until-culprit headline** (validated off-target core only).
- `agent.py`: pass the outcome panel into the evidence packet under a distinct
  `model_predicted_modules` key and instruct the narrative to mention liver/mito/reactive
  findings **only as lower-tier, model-predicted, needing confirmation** — still grounded, no
  fabrication. Keep the graceful templated fallback (extend `_fallback` to add one lower-confidence outcome sentence).
- **Acceptance (`tests/test_app.py`):** AppTest on troglitazone → the "model-predicted" section
  is present and mentions hepatotox; the assays-until-culprit metric/headline still renders and is
  off-target-derived (unchanged by outcome modules). Narrative for troglitazone (fallback path, no API key) contains
  "model-predicted" and does not claim a probability. Existing render tests stay green.

### [x] N7 — Mechanism-edge network viz (STRETCH — only if N1–N6 done + time)
- `render.py`: `mechanism_graph_dot(result, plan) -> str` — a Graphviz DOT string with the
  candidate at center and each flagged target's linked failed drugs as nodes, edges labelled by
  target/mechanism, node color by organ system (cardiac/CNS/GI/autonomic). Render with
  `st.graphviz_chart(dot)` (no local graphviz binary needed — frontend renders DOT). If
  `st.graphviz_chart` errors, fall back to a `networkx`+`matplotlib` PNG via `st.image` — but
  **do not add networkx/matplotlib unless the graphviz path fails in AppTest** (matplotlib is
  already present via streamlit? verify; if absent, keep graphviz-only and skip the fallback).
- **Acceptance:** AppTest on rimonabant → a graphviz chart element is present; the DOT string
  contains "rimonabant" or "taranabant" and "CB1".

---

## 6. HOW TO OPERATE THIS LOOP (each iteration)
1. Read this file; find the first `[ ]` milestone.
2. Implement it fully. Reuse existing code (`core._Engine`, `exp_score.loo_exclude`,
   `reference_failures.json` pattern). Do not rewrite the off-target class-matching math.
3. Run its acceptance check with `.venv/bin/python -m pytest app/tests -q` (add the new test
   file first) and/or the AppTest snippet. UI milestones: verify via `streamlit.testing.v1.AppTest`.
4. If it passes: flip `[ ]`→`[x]` here, commit (`N<n>: <what>`), continue.
5. If blocked (e.g. a curated SMILES won't parse): fix the data, or note it and implement the
   graceful path; never leave the suite red.
6. **Stop** when N1–N6 are `[x]` and `app/tests` is green. N7 is optional. Do not gold-plate.

## 7. DO-NOT list
- Do not merge outcome modules into the assays-until-culprit headline or present them above the
  validated off-target plan — they are tier-2.
- Do not present any score as a probability of harm.
- Do not let leave-one-out demo mode affect the default/novel path (default OFF; removes nothing when off).
- Do not add new pip dependencies (TDC/ProTox stay optional stretches behind graceful guards;
  none of N1–N7 require them).
- Do not invent drugs, mechanisms, citations, or validation numbers — cite FINDINGS.md / §4 / the
  curated-set sources; verify every curated SMILES parses.
- Do not commit `panel_actives.json` or other regenerable caches (already gitignored). The
  curated `outcome_actives.json` / `reactive_alerts.json` ARE committed (they are authored data,
  not fetched caches).
```
