# Toxicity Risk → Assay Recommender

**A tool that reorders a small-molecule candidate's in-vitro safety-assay plan to minimize *time- and cost-to-kill* — it moves the experiment most likely to catch a real, program-ending liability to the front of the queue, by linking the candidate to drugs that failed in the clinic for toxicity via a shared off-target mechanism. Validated on 20 historical failures: for off-target liabilities a standard panel buries, it cuts the number of assays you run before hitting the culprit from ~11 to ~4.**

> Positioning in one line: *"For a novel candidate that is going to die on a safety liability, we surface the assay that will kill it **first** — so you spend 4 experiments reaching the kill decision instead of 11. We do it by connecting the candidate to clinically-failed drugs it shares a toxicity **mechanism** with — including ones no structure search would find. It prioritizes experiments to fail cheaper and faster; it does not declare a molecule safe."*

---

## 0. Project context & status (read first)

- **Event:** a bio × AI hackathon. Deliverable is an **MVP + investor/judge-facing demo** (§6), not a validated product.
- **Core method:** **off-target class matching**. For each safety-panel target, score the candidate against the whole ChEMBL active class for that target, then z-score it against a background set. This replaced the early research shorthand and should be the term used in product/research docs.
- **The core value thesis — fail faster, fail cheaper.** In drug discovery most candidates die, and the expensive waste is *how long and how much you spend before the kill decision*. A broad in-vitro secondary-pharmacology panel is run largely in parallel or in a generic default order. If a candidate is doomed by a specific off-target liability, the assay that proves it is often **not** the one you run first. This tool reorders the plan so the program-ending experiment comes first → **less time, less spend, fewer animals, earlier go/no-go.** The product metric is **assays until culprit** (how many experiments until you hit the real liability), and we lower it.
- **De-risking is DONE — 8 experiments, all on real ChEMBL data** (scripts + results in **[`experiments/`](experiments/)** and **[`experiments/derisk/`](experiments/derisk/)**). Headline results below; details in §4. The verdict: **build it — the value is real, but narrower and more honest than a first pass would claim.**
- **What the de-risking changed about the pitch (all forced by our own data):**
  1. **The value is *buried* off-target liabilities, explicitly NOT hERG.** A standard panel already front-loads hERG; there we add nothing and are sometimes worse. Our measured value is on liabilities a default panel buries at rank 11–15 — **5-HT2B valvulopathy, CB1 psychiatric tox, 5-HT3 colitis** — which we lift to rank 1. Do **not** lead with the terfenadine/hERG story; keep it only as a clean recovery-validation anchor.
  2. **The engine's advantage is conditional on *novelty*.** Against a candidate that already resembles a known failure, cheap 2D similarity matches us. Our engine earns its keep — and wins decisively — on **novel chemotypes with no analog in the failure database**, which is precisely the real use case. State the claim conditionally; a judge who runs the obvious ablation will otherwise catch an overclaim.
  3. **The addressable population is *narrow* — the value metric does not generalize to all withdrawn drugs.** A blind hold-out (drugs *and* culprit labels selected by ChEMBL, not us — §4) confirmed the *magnitude* (assays until culprit 4.7 vs 11.4) but exposed the scope: when an objective source defines the culprit, most withdrawn drugs' dominant panel activity is **on-target pharmacology** or their withdrawal cause is **off-panel** (hepatotoxicity dominates: 224/1,014 warnings, a blind spot for off-target class matching). Drugs withdrawn *specifically for a reachable secondary-pharmacology off-target* (the pergolide/rimonabant/fen-phen class) are a **real but narrow slice**. Frame the value as *"for candidates whose failure mechanism is a secondary-pharmacology off-target"* — **not** "assays until culprit 11→4 across withdrawn drugs" broadly.
- **Repo state:** the working demo is in **[`app/`](app/)**, with tests in **[`app/tests/`](app/tests/)**. This doc is the long research memo; **[`SUMMARY.md`](SUMMARY.md)** is the one-page version; **[`experiments/`](experiments/)** contains the feasibility scripts and the 8 de-risking tracks under `experiments/derisk/`, each with a `FINDINGS.md`.
- **Parameters still TBD:** timeframe/deadline; team size & strengths; explicit judging criteria; target LLM + infra.
- **Data caveat:** dataset sizes are **last-published figures** and several conflict across mirrors (e.g. ClinTox is **1,484** on TDC and **1,491** in the MoleculeNet paper). **Re-confirm every count and DOI against the primary source before quoting in a pitch.**

---

## 1. The idea

- **Input:** a small-molecule candidate (SMILES; optionally a 3D conformer).
- **Reference database:** drugs annotated with **why they failed** — prioritizing those that **failed clinical trials or were withdrawn specifically for a toxicity**, each tagged with the **organ, phenotype, and mechanism** and a provenance link.
- **Output — a reordered in-vitro assay plan optimized for time-and-cost-to-kill:**
  1. **The headline — a re-ordered safety-assay plan:** the assay most likely to catch a real, program-ending liability moved to the front, each tagged with an **actionability label** (*no-go / counter-screen / monitor*, §4 severity/actionability weighting) and an **evidence tier** (*mechanism-linked / lower-confidence model prediction / abstain*). The value is quantified as **assays until culprit** — how far down the default plan the killer experiment sits vs. where we put it. E.g. *"5-HT2B counter-screen first (no-go: resembles drugs withdrawn for valvulopathy) — default panels bury it at rank 11; we move it to rank 1."*
  2. A **per-liability priority profile** (supporting) — enrichment vs. a matched lower-liability comparator, or a 0–100 priority score; **not** a probability of harm (§4 Layer 3).
  3. An **evidence trail** — the linked failed drugs, the organ/mechanism, evidence tier, and provenance.

Under the hood it's a **tiered portfolio** (§4): the validated off-target core plus lower-confidence outcome modules for the metabolism-driven liver/mito gap. It is **decision-support / experiment-prioritization**, not a binary safety classifier — more honest (proxies miss metabolism, dose, exposure) and directly aligned with the economic goal: **reach the kill decision with fewer experiments.**

**Who it is for, and when it helps.** The tool is for a team facing a broad, expensive safety panel on a **novel** candidate. It helps most exactly when the candidate does *not* obviously resemble a known bad drug (if it did, you'd catch the liability by eye and wouldn't need us). It is a triage layer on top of, not a replacement for, the assays themselves.

---

## 2. Why it's novel (white space + prior art to differentiate from)

The three halves each exist separately; the **combination is unoccupied** (documented, reproducible search in **[`novelty-search-log.md`](novelty-search-log.md)** — 11 structured queries across Amass/PubMed/web, 2026-07-04):

- **Endpoint predictors** — ProTox 3.0, admetSAR 3.0 / ADMETlab 3.0, ADMET-AI, DeepTox, DEREK. Predict *tox probabilities / structural alerts*, **not "which assay to run,"** and are **not grounded in clinical-trial failures**. (We measured ProTox directly — see §4 competitor check.)
- **Read-across engines** — EPA **GenRA**, OECD QSAR Toolbox, AMBIT. Analog-based tox inference, but reference sets are **environmental chemicals** and output isn't a prioritized clinical assay panel.
- **Off-target liability engines** — **AbbVie OTSA** (in-silico secondary-pharmacology over >7,000 targets, ~2M preclinical+clinical alerts) is the closest *engine* — essentially off-target class matching plus an alert layer. It predicts liabilities but does **not** output a ranked in-vitro assay plan grounded in a curated clinical-failure set. **ToxEvaluator** (Pfizer 2016) is the closest A×B integration (similarity + off-target + tox linkage) but is post-hoc interpretation, not a prospective ranked recommender.
- **Closest paper** — Albert, Skinner et al., *"Comprehensive Analysis of Clinically Discontinued Compounds Using an In Vitro Secondary Pharmacology Panel…"* (**ACS Pharmacol Transl Sci 2025, DOI 10.1021/acsptsci.5c00452; PMID 40969870** — verified): **52 compounds discontinued 2001–2021 for safety** run through the **SAFETYscan47** (47-target) panel. Confirmed: it screens **everything** and argues the primary target doesn't predict off-target risk — it does **not** use similarity to *prioritize a subset* of assays. **That gap is our niche.**

**Our defensible novelty — an integration, not a new capability.** Every *pairwise* edge is occupied; two items (OTSA, Albert 2025) sit one step away. The white space that survived the search is the **specific deliverable**: a **ranked, clinical-failure-grounded, severity-weighted assay plan** with an **assays until culprit** value metric — that exact output appeared in **zero** sources. Lead with the deliverable, not with "we invented off-target prediction."
1. **Reference set grounded in *toxicity-causal* clinical failures** (not assay-positives, not mere market withdrawals).
2. **Output = a specific, ranked assay panel** optimized for cost-to-kill, made credible by organ+mechanism-resolved labels and the validation in §4.

> **Be precise about label strength in the pitch:** the *clinical-failure grounding* is what's novel, but the *scoring backbone* also uses organ-liability labels from DILIrank/DICTrank plus off-target class matching. Don't imply the whole score comes from toxicity-causal clinical failures: curated failure labels provide the clinical grounding, while organ-liability labels and off-target class matching do much of the statistical lifting.

> **Must-differentiate in the deck** (a judge may know these): Albert 2025 (SAFETYscan47), AbbVie OTSA, ToxEvaluator (Pfizer 2016), Deaton 2019 (panel prioritization by genetics), Zatorski 2025 / Xue 2013 (withdrawn-drug similarity → tox flag).

---

## 3. Datasets

### 3a. Toxicity-causal failure set (the novelty; assemble + curate, don't just merge)
No public dataset cleanly says "these N molecules died in the clinic *because of toxicity*, here are their SMILES + organ + mechanism." Build it — and **curate the *reason*, don't assume it**:

| Dataset | Content | Size | Access | Role |
|---|---|---|---|---|
| **ClinTox** (MoleculeNet / TDC) | FDA-approved vs. **failed clinical trials for toxicity** | **1,484 (TDC) / 1,491 (paper)**; tox-failure positives only ~90–110 | Free — `deepchem`, TDC, HF `zpn/clintox` | Start here. Positives mined from AACT. Class-imbalanced. |
| **WITHDRAWN 2.0** (Charité) | Withdrawn/discontinued drugs + withdrawal reasons | ~578 (v1) → ~631–635 (2.0) | Free — cheminfo.charite.de/withdrawn; NAR 2024 52:D1503 | ⚠️ **Includes *non-safety* withdrawals**. Must filter on `safety_related` + organ before calling anything a curated toxicity-causal failure. |
| **Probes & Drugs** — withdrawn set | Withdrawn drugs pre-joined to ChEMBL/DrugBank | ~633 (verify release/date) | Free — probes-drugs.org/compoundsets | De-dupe / canonical-SMILES join. |
| **DrugBank** (withdrawn flag) | `withdrawn` category, structures | ~250–300 (verify) | Free academic; **non-commercial** | Structure + flag; weak on reason. |
| **AACT** (ClinicalTrials.gov dump) | `why_stopped` free-text, `TERMINATED/…`, interventions | full CT.gov | Free — aact.ctti-clinicaltrials.org | **Stretch.** NLP-classify tox terminations → resolve names to SMILES (~30–60%). Demo a slice; don't rebuild wholesale. |

> **Note on curation:** structures resolve cheaply (PubChem name→SMILES: 20/20 in our expanded set), but the *reason / organ / mechanism / tier* fields are a manual sink (~2 person-weeks for ~200 drugs; see §4 Phase-0 curation). For the MVP the reference set is small and hand-curated (the drugs in `experiments/derisk/drugs.json` + fixtures).

### 3b. Organ + mechanism-resolved toxicity labels (power the "which assay" step)
FDA-curated *ranked* lists — highest-quality organ-resolved labels; backbone of scoring **and** of the matched comparator (§3e).

| Dataset | Content | Size | Assay it points to |
|---|---|---|---|
| **DILIrank** (FDA LTKB) | **liver injury**, 4 classes incl. **No-DILI-concern** | ~1,036 (v1); **2.0 = 1,336** | Hepatocyte/spheroid, BSEP inhibition, mito-tox |
| **DILIst** (FDA) | binary hepatotoxic (broader) | ~1,279 | same |
| **DICTrank** (FDA, 2023) | **cardiotoxicity**, 4 concern levels incl. **No-DICT-concern** | ~1,318 | hERG, iPSC-cardiomyocyte, Nav1.5/Cav1.2 · DOI 10.1016/j.drudis.2023.103770 |
| **hERG datasets** | hERG channel block | TDC 648; DeepHIT ~14k; UnihERG ~20k | hERG patch-clamp / electrophysiology |
| **SIDER** | side effects from labels | 1,556 drugs | phenotype → organ → assay |
| **OFFSIDES / TWOSIDES** | FAERS-derived side effects / DDIs | 438,801 / 868k assocs | pharmacovigilance layer |
| **LiverTox** (NIH) | narrative hepatotox monographs | ~1,200 agents | qualitative liver context (text) |

> The **No-concern** classes in DILIrank/DICTrank are the best public **lower-liability comparators** (§3e) — better than "any approved drug."

### 3c. Assay / mechanism datasets (map mechanism → a named assay)
| Dataset | Content | Size | Role |
|---|---|---|---|
| **Tox21** | 12 pathway assays | ~7,800 × 12 | each pathway *is* an assay |
| **ToxCast** (EPA) | HTS in-vitro endpoints | ~800–1,400 endpoints × >4,200 chem | richest compound → assay-endpoint matrix |
| **ChEMBL** | measured bioactivities incl. safety targets (hERG, 5-HT2B, BSEP…) | ~2.4M cpds | **the mechanism layer** — reference drugs' *known* off-targets; MCP available; **off-target class matching runs on this** |
| TG-GATEs / DrugMatrix / LINCS | toxicogenomics / histopath / cell response | ~170 / ~600 / >20k | mechanism enrichment (stretch) |

### 3d. Structure sources
**PubChem** (name→structure resolution — verified reliable, note API now returns `SMILES`/`ConnectivitySMILES` properties), **ChEMBL** (CC-BY-SA; structures + bioactivity; MCP), **DrugBank** (annotations; non-commercial).

### 3e. Matched lower-liability comparator (the denominator — do not call it "safe")
The class-conditional score needs a denominator. **"Approved" ≠ "safe."** Use **per-organ No-concern classes** (DILIrank *No-DILI-concern*, DICTrank *No-DICT-concern*) as the comparator, matched to the toxic set by scaffold / physchem / route where possible. In the de-risking we used a fixed 24-drug background set (see `experiments/score.py`) as a lightweight stand-in; a production build should upgrade to the matched No-concern comparator.

---

## 4. Architecture

### The engine's job, stated as the value metric
For a candidate, rank the safety-panel assays by predicted-liability so the **program-ending assay is first**. Success = **low assays until culprit** (rank of the killer assay) on held-out real failures, vs. a defensible default panel order. Everything below serves that metric.

### Empirical reality check #1 — 2D structure can't link mechanistic failures (drove the design)
7 pairs of drugs **both withdrawn/failed for the same organ toxicity via the same documented mechanism**; 2D fingerprints miss all of them (ECFP4 0.05–0.19; feature/pharmacophore FCFP only 1.3–1.5× lift). The signal linking them is a **shared binding target**, invisible to 2D structure. USRCAT 3D shape was tested and is a **weak supplement** (helped 4/7, hurt 3/7 — global shape dominates over the target-relevant toxicophore). **Conclusion: the retrieval engine must match off-target mechanisms, not just 2D structure.** (Scripts: `run_all.py`, `run_pharm.py`, `usrcat_experiment.py`.)

### Empirical reality check #2 — off-target class matching recovers the mechanism (the fix)
For each safety target, fetch its known actives from ChEMBL and score a query by mean top-5 Tanimoto to that *class* (strict leave-one-out: query + named partners removed first). **terfenadine → hERG = 0.57, z = +6.4** vs pairwise terfenadine–cisapride = 0.19 — and cisapride wasn't even in the fetched hERG set, so recovery is driven entirely by *other* hERG ligands. Across 6 query drugs the correct mechanism was the top-scoring target for 4/6 and above background for 6/6. **Why it works:** terfenadine isn't similar to cisapride, but each is similar to *some* member of the large, diverse set of known hERG ligands — aggregating over the whole target class is far more sensitive than one pairwise comparison. (Scripts: `fetch_actives.py`, `score.py`.)

### The validation that decided go/no-go — assay-recovery at scale (n=20)
*(`experiments/derisk/expanded/` — the decisive experiment. Full data in its `FINDINGS.md`.)*

20 withdrawn/failed drugs, each scored against an **18-target safety panel**, strict leave-one-out (query + all same-culprit mechanistic partners removed from every ChEMBL class first). Killer-assay rank recorded under **Ours** (rank by candidate z-score) vs **Baseline B** (defensible, candidate-independent standard-practice order: cardiac ion channels front-loaded, buried off-targets low) vs **Baseline C** (data-derived base-rate order). Split into **buried** (culprit = a non-cardiac-channel off-target — the real test, n=10) and **cardiac_hERG** (contrast, n=10).

| | Top-3 recovery (Ours / B / C) | Mean **assays until culprit** (Ours / B) |
|---|---|---|
| **Buried (n=10)** | **7 / 1 / 5** | **3.8 / 11.3** |
| **Cardiac_hERG (n=10)** | 5 / 10 / 10 | 3.6 / 1.0 |

**★ The value number:** on buried off-target liabilities, we move the killer assay from a mean rank of **11.3 → 3.8** — roughly a **3× reduction in assays until culprit**, i.e. in experiments (time, cost, animals) spent before the kill decision. **6/10 are genuine non-obvious wins** — Ours puts the killer in top-3 *and* the default panel does not — spanning **three distinct mechanisms**: 5-HT2B valvulopathy (pergolide, cabergoline, methysergide), CB1 psychiatric tox (rimonabant z=+10.8, taranabant z=+9.6), 5-HT3 colitis (alosetron). The earlier n=6 result (one win, pergolide) generalized.

> **Scope qualifier (from the blind hold-out below — read before quoting the number).** This 11.3→3.8 is measured on cases where the culprit *is* a secondary-pharmacology off-target we can reach. It is **not** a claim about withdrawn drugs in general — most are withdrawn for off-panel liver/hematologic/idiosyncratic tox, or their dominant panel activity is on-target pharmacology. Quote the number as *"for candidates whose failure mechanism is an off-target liability,"* never unqualified.

### Blind hold-out — the number survives external selection, but the addressable class is narrow
*(`experiments/derisk/blind_holdout/` — the rigor check against the "you cherry-picked" critique.)*

To remove human selection bias, **ChEMBL selected the drugs** (withdrawn `drug_warning` table → 224 safety withdrawals → 20 with a measured panel culprit, forced disjoint from our hand-picked set) **and ChEMBL labeled the culprit** (most-potent measured pChEMBL≥6 panel target). Result: buried top-3 **12/19 vs default 0/19**, assays until culprit **4.68 vs 11.37** — the magnitude **replicates** on a set we neither picked nor labeled, so the cherry-pick critique on *selection* is defeated.

**But external labeling revealed the honest boundary:** only **1/20 drugs** is cleanly anchored to its actual withdrawal cause; for the other 17 the "objective culprit" is the drug's **on-target pharmacology** (4 NSAIDs→COX-1, remoxipride→D2, indalpine→SERT) or its true cause is **off-panel** (nefazodone/troglitazone/ketoconazole → withdrawn for hepatotoxicity, which off-target class matching cannot reach). The one clean buried off-target surfaced on an unpicked drug: **lorcaserin → 5-HT2B**. **Takeaway:** the method is valid and the default panel genuinely buries non-cardiac liabilities, but the class of "withdrawn *for* a reachable off-target" is a narrow slice of all withdrawals — position accordingly.

**The honest boundary, measured not asserted:**
- **On hERG/cardiac we add nothing and are sometimes worse** (Ours 5/10 vs default 10/10). Real panels already front-load hERG at rank 1; for phenothiazines the drug's own D2 on-target outranks hERG, and for fluoroquinolones ligand-similarity doesn't capture the atypical hERG block. **The tool's value is explicitly not hERG.**
- **The only misses are metabolite-active drugs** — fenfluramine, dexfenfluramine (rank 13), partial sibutramine — where the toxic species is a metabolite the parent structure can't reveal. The documented boundary held exactly.

### Is off-target class matching even necessary? — the ablation (adversarial)
*(`experiments/derisk/ablation/` — off-target class matching vs a dumb "2D nearest-failed-drug" baseline on the same 20-drug assay-ranking task.)*

**Taken at face value, naive 2D wins overall (16/20 vs 12/20 top-3).** This is the most important nuance in the whole project, and it sharpens the claim:
- Naive's win is largely **artifact**: hERG is half the reference set (so "nearest failed drug" defaults to hERG for any drug-like query), and several wins ride on **near-duplicate twins** (fenfluramine↔dexfenfluramine sim 1.00; thioridazine↔mesoridazine 0.70; the two fluoroquinolones 0.50) — leakage-adjacent, not generalization.
- **Under strict leave-one-out (simulating a genuinely novel candidate with no analog in the failure DB — the actual use case), naive 2D collapses to 0/20 while off-target class matching holds 12/20.** If a candidate already resembles a known failure you don't need *any* tool; the tool exists for the novel molecule, and that is exactly the regime where only class aggregation works.
- Off-target class matching **uniquely wins the singleton/diverse mechanisms** — taranabant (only 0.15 similar to its lone CB1 partner, yet class matching nails it at z=8.6, rank 1), alosetron (5-HT3), mibefradil (Cav1.2). An aggregation control (compare to *one* random ChEMBL active instead of the class) scores just 2/20 — confirming it's the **diverse-class aggregation** that carries the signal, not merely "having ChEMBL data."

**The claim this licenses (state it exactly, don't overreach):** *"For a novel candidate with no analog among known failures, pairwise similarity buries the liability at rank ~12 and only off-target class matching recovers it to the top (naive 0/20 → class matching 12/20 under strict hold-out)."* Do **not** claim class matching beats similarity search in general — on a reference set with analogs present, it doesn't. This makes the applicability-domain check (below) load-bearing: the tool should **route** — *known analog → cheap similarity suffices; novel → class matching adds the value.*

### Competitor check — "why not just use ProTox?"
*(`experiments/derisk/protox/` — a real REST API was reverse-engineered and a runnable driver written; live run was IP-throttled, so the verdict is analytical from ProTox-3.0's full endpoint inventory.)*

ProTox-3.0 emits **organ/pathway/target *probabilities*** (13 organ models incl. a hERG-trained `cardiotoxicity` flag, 12 Tox21 pathways, 14 MIEs, 6 CYPs, and 15 "toxicity targets" gated on a **solved human crystal structure** — 15 of 73 panel targets). Verdict:
- **Cardiac:** its cardiotoxicity model fires — ~matches the *easy* cases, as a binary organ flag (not "run hERG first").
- **Buried:** **structurally cannot cover 5-HT3, NET, or Cav1.2** (no human crystal structures); 5-HT2B/CB1 *might* be among its 15 (unconfirmed). Even then it's a flat pharmacophore-fit score.
- **Never** produces a ranked, reordered assay panel, and **never** grounds a signal in the specific clinically-failed drugs sharing the mechanism.

**Blunt answer:** ProTox **partially overlaps** on coarse cardiac tox and is a legitimate lower-confidence outcome component (we already wrap it) — but it is **not a substitute** for the linkage-and-assay-triage layer, and it does not make the project redundant. Its structure-gated 15-target ceiling is exactly why a ligand-similarity engine (needs no structure; reaches 37/39 panel targets) is the right tool for the buried liabilities.

### De-risking summary — the four earlier cheap checks (all run; `experiments/`)
| Check | Result | Verdict |
|---|---|---|
| **Specificity** (hERG) | ROC-AUC **0.894**, blockers vs assay-confirmed non-blockers; 19/19 common drugs low | ✅ discriminates; no false-positive flood |
| **Cross-target replication** | AUC **SERT 0.95, AChE 0.92, MAO-A 0.88** (hERG 0.89 replicates) | ✅ holds across GPCR/transporter/enzyme classes |
| **Scaffold generalization** | Pooled leave-one-scaffold-out AUC **0.913**; diverse single-ring holdout **0.667** | ⚠️ real, but degrades for truly novel isolated chemotypes |
| **Panel coverage** | **37/39 targets serviceable**; gap: **BSEP 5 actives; mito = no single target** | ✅ broad off-target reach; ❌ liver is a data desert → lower-confidence outcome modules |

### The tool is a *portfolio*, not one engine
Breadth comes from the right method per domain, fused into one prioritized assay panel, each item tagged with an **evidence tier**:

| Module | Domains | Method | Evidence tier (shown in UI) |
|---|---|---|---|
| **Validated off-target core** | cardiac, CNS, GI, endocrine, immune (37 targets) | similarity to known ChEMBL actives → linked clinical failures | **mechanism-linked** (strongest) |
| **Lower-confidence outcome modules** | liver DILI, mitochondrial, reactive-metabolite | QSAR/read-across on DILIrank, Tox21-mito, structural alerts; wrap ADMET-AI / ProTox | **model-predicted / alert-based** |
| **Fusion** | all | rank + tag every recommendation by evidence tier & confidence; **abstain** or flag a *known analog* outside domain | — |

The novelty is the **fusion + clinical-failure grounding + assay triage + honest evidence-tiering + the cost-to-kill framing**, not any single predictor. Coverage flags are first-class output: *covered mechanism / weak coverage / abstain / known analog*.

### Pipeline
```
             ┌─────────────────────────────────────────────┐
  SMILES  →  │ 1. Multi-method similarity / linkage         │
             │   2D fingerprint screen (ECFP, FPSim2)       │  ← fast pre-filter, obvious analogs
             │   feature/pharmacophore fingerprint (FCFP)   │
             │   3D shape + pharmacophore (USRCAT)          │  ← test-gated, weak supplement
             │   off-target class matching (ChEMBL actives) │  ← differentiator for novel chemotypes
             └───────────────────┬─────────────────────────┘
                                 ▼
             ┌─────────────────────────────────────────────┐
             │ 2. Link to failed-drug reference DB          │
             │   • candidate ⇄ drugs w/ shared mechanism    │
             │   • carry organ + phenotype + tier + source  │
             └───────────────────┬─────────────────────────┘
                                 ▼
             ┌─────────────────────────────────────────────┐
             │ 3. Per-liability scoring (class-conditional) │
             │   • enrichment vs matched No-concern class   │
             │   • z-score vs background; N× or 0–100       │
             │   • applicability-domain / known analog flag │
             └───────────────────┬─────────────────────────┘
                                 ▼
             ┌─────────────────────────────────────────────┐
             │ 3b. Severity / actionability re-weight       │
             │   • no-go / counter-screen / monitor         │
             └───────────────────┬─────────────────────────┘
                                 ▼
             ┌─────────────────────────────────────────────┐
             │ 4. mechanism → assay mapping → REORDERED plan│
             │   → assays ranked to minimize assays-to-kill │
             │   → lower-confidence liver/mito modules fuse │
             └───────────────────┬─────────────────────────┘
                                 ▼
             ┌─────────────────────────────────────────────┐
             │ 5. Explainability + Agent + Visualization    │
             └─────────────────────────────────────────────┘
```

### Reference DB schema (split labels — avoid circular recommendations)
Keep organ, phenotype, mechanism, and assay as **separate fields** so the recommended assay is *derived*, not baked into the label:
```
compound_id | canonical_SMILES | source(s) | provenance/evidence_link
failure_reason (free text) | safety_related (bool) | fail_stage | label_tier
organ_system | tox_phenotype | mechanism (e.g. hERG block, BSEP inhib, 5-HT2B agonism)
known_offtargets (from ChEMBL) | assay_endpoint (observed) | recommended_assay (derived)
```
- **`safety_related` + `failure_reason` are mandatory** — never treat a compound as a curated toxicity-causal failure without a reason + provenance.
- **Label strength:** curated toxicity-causal failures are strongest; organ-severity-ranked marketed-drug labels from DILIrank/DICTrank are useful but weaker; broad side-effect or label-derived flags are noisiest.

### Layer 1 — off-target class matching (the safety-panel pipeline)
The mechanism vocabulary is a **defined secondary-pharmacology safety panel**: the **Bowes et al. ~44-target minimal panel** (*Nat Rev Drug Discov* 2012;11:909), aligned to **SAFETYscan47** (Albert 2025). Each off-target → an **organ** and a **named assay** (Layer 4). Pipeline:
- **(a) Reference side** — annotate each clinically-failed drug with its **measured** ChEMBL panel activity + organ/mechanism/provenance. Facts, not predictions.
- **(b) Candidate side** — score against the **same panel**: for each target, mean top-5 Tanimoto to its ChEMBL actives (pChEMBL ≥ 6), as a **z-score vs a background set of non-binder drugs** (prototype: terfenadine→hERG z = +6.4). Production upgrade: swap raw class-Tanimoto for a trained target-prediction model (SEA / per-target QSAR).
- **(c) Link** — candidate's predicted panel-hits → failed drugs sharing them → organ → assay → reordered plan.
- **Fallbacks:** target with too few actives (BSEP = 5) → fall back to direct reference-drug annotations or the 2D floor, and flag low coverage. **Metabolite-active drugs:** score parent AND known active metabolite(s), aggregate by **MAX** (§ metabolite rule below).

### Layer 3 — Per-liability scoring (class-conditional; enrichment, not a probability)
Similarity to a toxic drug means nothing unless the candidate is *more* linked to toxic drugs than to the **matched No-concern comparator** (§3e). Report as **"4.1× liver-liability enrichment"** or a **0–100 priority score** — never `0.71`-style decimals that read as P(harm) unless calibrated. Attach: driving links (drug, organ, mechanism, tier, provenance) + **applicability-domain / known analog flag**.

### Layer 3b — Severity / actionability weighting (likelihood ≠ consequence)
The score measures *how likely* the candidate hits a target — **not how much that hit matters for the kill decision.** Ranking must weight **consequence**:
- **Per-target severity tier** (hand-authored, Bowes-2012-derived — starter table below): hERG / 5-HT2B = high (potential no-go); H1 = low (monitor).
- **Clinical-failure grounding *is* a severity signal** — resembling a curated clinical failure at a target ("this hit killed real programs at clinical exposure") outranks resembling a benign ligand. The failure label carries both mechanism and consequence.
- **Exposure margin** — off-target potency vs efficacious exposure; for a novel candidate absolute exposure is unknown → flag *"margin needs exposure context,"* ground via reference failures' known margins where available.

Ranking ≈ (hit confidence) × (target severity) × (clinical-failure grounding) × (margin, if known). Tag every recommendation **no-go / counter-screen / monitor**.

**Starter severity/actionability table** (representative — expert-review before the pitch):

| Off-target | Organ / phenotype | Severity | Default action |
|---|---|---|---|
| hERG (KCNH2) | cardiac / QT, torsades | **high** | early-kill flag → hERG patch-clamp + iPSC-CM *(note: default panels already front-load this — low marginal value)* |
| **5-HT2B** | cardiac valve / valvulopathy | **high** | **early-kill flag → 5-HT2B counter-screen** *(buried by default panels — high marginal value)* |
| **CB1** | CNS / psychiatric, suicidality | **high** | **early-kill flag → CB1 counter-screen** *(buried — high marginal value)* |
| **5-HT3** | GI / ischemic colitis | high | **5-HT3 counter-screen** *(buried — high marginal value)* |
| Nav1.5 / Cav1.2 | cardiac / conduction, contractility | high | counter-screen |
| µ-opioid | CNS / respiratory depression, abuse | high | counter-screen + exposure margin |
| MAO-A | CV/CNS / hypertensive crisis | high | counter-screen |
| AChE | autonomic / cholinergic crisis | high | counter-screen |
| PDE3 | cardiac / inotropy, ↑ mortality | high | counter-screen |
| D2 | CNS / EPS, hyperprolactinemia | med | counter-screen |
| M1–M3 muscarinic | autonomic / anticholinergic | med | monitor / counter-screen |
| α1A adrenergic | CV / orthostatic hypotension | med | monitor |
| COX-1 | GI / bleeding | med | monitor |
| SERT / NET | CNS/CV / serotonin syndrome, BP | med | monitor |
| H1 | CNS / sedation | **low** | monitor |

### Layer 4 — mechanism → assay mapping (derived, hand-authored; carries the Layer-3b tag)
| Mechanism / organ | Recommended assays (priority order) |
|---|---|
| **5-HT2B agonism / valvulopathy** | **5-HT2B binding/functional counter-screen** (per Setola 2003) |
| **CB1 antagonism / psychiatric** | **CB1 binding/functional counter-screen** + CNS safety pharmacology |
| **5-HT3 / ischemic colitis** | **5-HT3 counter-screen** + GI motility |
| hERG block / cardiac (DICTrank) | hERG patch-clamp → iPSC-cardiomyocyte → Nav1.5/Cav1.2 |
| BSEP inhibition / cholestatic liver | BSEP inhibition → hepatocyte spheroid |
| Mitochondrial / hepatocellular liver | mito assays (Glu/Gal, Seahorse) → hepatocyte spheroid |
| Reactive metabolite | GSH-trapping / covalent-binding; CYP-dependent cytotox |
| Tox21/ToxCast pathway hit | the specific named HTS assay for that pathway |

### Layer 5 — Explainability, agent, visualization
- **Attribution:** name the *shared mechanism/off-target*, the linked failed drugs, and where structural, the driving substructure (fingerprint-bit→atom mapping for 2D fingerprints; shared pharmacophore for pharmacophore/shape methods). Defends against "isn't this just fingerprint similarity?"
- **LLM agent** turns retrieved evidence into a report (retrieval + reasoning; no training for v1).
- **Visualization** — candidate at center, linked failed drugs around it, **edges labeled by mechanism**, colored by organ; the reordered assay plan with the **assays until culprit** delta shown explicitly (default rank → our rank).

### Metabolite rule (from de-risking)
*(`experiments/derisk/metabolite/`)* Score the parent AND known active metabolite(s); aggregate by **MAX** (never mean). MAX lets a toxic metabolite rescue an inactive parent (fenfluramine +0.35→+1.22; benfluorex +0.64→+1.22) while a detoxified metabolite of a toxic parent cannot dilute the parent's signal (terfenadine→fexofenadine stays +6.38). **Limit:** partial rescue only (still a weak positive), and it requires external metabolism knowledge unavailable for a novel candidate without CYP site-of-metabolism prediction.

### Applicability domain / abstain / known analog (from de-risking)
*(`experiments/derisk/abstain/`)* The intuitive nearest-neighbour rule **fails** (in-domain drugs and OOD molecules overlap in max-NN Tanimoto). The rule that works is a **descriptor-box scope gate**: abstain if it doesn't parse / <6 heavy atoms / contains a metal, OR MW > 800 / fraction-halogen > 0.30 / fraction-oxygen > 0.33. Measured **0/31 false-abstain, 0/11 false-accept** on the test set. Report a low best-z as `weak_coverage`, distinct from `abstain`. **Add a known analog flag** (candidate has a high-similarity match among known failures) → tells the user cheap similarity already suffices and off-target class matching adds little (per the ablation). Demo an abstain with a **wide-margin** molecule (cyclosporine, MW 1203); avoid borderline cases.

### Rigor: validation, calibration & applicability domain
- **★ Assays until culprit — THE headline metric (measured, n=20):** buried liabilities **11.3 → 3.8**; top-3 recovery **7/10 vs 1/10** default. This *is* the product promise (reach the kill decision with fewer experiments). Target-AUC and organ discrimination are **supporting**.
- **Scaffold-split CV (never random)** — measured pooled AUC 0.913, novel-isolated 0.667.
- **Per-target discrimination (supporting):** hERG 0.89, SERT 0.95, AChE 0.92, MAO-A 0.88 + enrichment factor / precision@k.
- **Necessity ablation (done):** off-target class matching vs naive 2D — class matching wins only in the novel-chemotype regime (strict leave-one-out naive 0/20 vs class matching 12/20); report it, don't hide it.
- **Baselines to beat:** naive 2D nearest-failed-drug (done), plain ECFP, and ≥1 published predictor (ProTox — coverage-analyzed; run the driver when un-throttled).
- **Calibration:** does 4× enrichment map to a real elevated failure rate? Report a calibration curve; ranked-but-uncalibrated is fine for prioritization if not sold as P(harm).
- **Honest bounds:** proxies miss metabolism/reactive-metabolites/dose; no value on hERG; off-target class matching only earns its keep on novel chemotypes; curated toxicity-causal failure labels are scarce.

---

## 5. Roadmap (revised after de-risking)

Phase 0 de-risking is **done and decisive** (8 experiments). The MVP and phase-2 demo are implemented in [`app/`](app/): the validated off-target core is framed as a cost-to-kill reorderer, and lower-confidence outcome modules add liver/mito/reactive-metabolite breadth. Remaining work is hardening, calibration, and scale-up.

**Phase 0 — De-risking (DONE ✅ — §4 + [`experiments/`](experiments/) + [`experiments/derisk/`](experiments/derisk/))**
- [x] 2D fails on mechanistic pairs; 3D/USRCAT a weak supplement.
- [x] Off-target class matching recovers the mechanism (terfenadine→hERG z=+6.4).
- [x] Specificity AUC 0.894; cross-target SERT 0.95 / AChE 0.92 / MAO-A 0.88; scaffold 0.913 pooled / 0.667 novel-isolated; coverage 37/39.
- [x] **Assay-recovery generalizes (n=20):** buried assays until culprit **11.3→3.8**, top-3 **7/10 vs 1/10**, **6 non-obvious wins** across 5-HT2B / CB1 / 5-HT3.
- [x] **Necessity ablation:** off-target class matching wins the novel-chemotype regime (naive 0/20 vs class matching 12/20 strict leave-one-out); ties/loses where an analog exists.
- [x] **ProTox competitor check:** partial cardiac overlap, structurally can't localize buried targets, no ranked panel → not redundant.
- [x] **Novelty search log:** triple-combination white space survives; differentiate OTSA / Albert 2025 / ToxEvaluator.
- [x] **Blind hold-out** (ChEMBL selects drugs *and* culprit labels): magnitude replicates (4.7 vs 11.4) but the addressable class is narrow — most withdrawals are off-panel/on-target. → qualify the value claim.
- [x] **Abstain rule:** descriptor-box gate works (0/31, 0/11); NN rule fails. **Metabolite rule:** MAX(parent, metabolite), partial rescue.
- → **Verdict: build it. Value is real, on buried off-target liabilities, for novel candidates. Honest about hERG (no value) and metabolite tox (blind).**

**Phase 1 — validated off-target core (IMPLEMENTED in [`app/`](app/))**
- [x] Reference set: SMILES + ChEMBL off-targets + curated failure fixtures in `app/data/reference_failures.json`.
- [x] Off-target class-matching module across the 18-target panel (class-membership z-scoring, live path without leave-one-out, demo/validation leave-one-out mode, MAX metabolite aggregation).
- [x] **Evidence-tier + coverage + known analog flags as first-class output.**
- [x] **Severity/actionability table** (§4 Layer 3b) with marginal-value emphasis (buried > hERG).
- [x] mechanism → assay mapping → **reordered plan with the assays-until-culprit delta** displayed.
- [x] **Assays-until-culprit metric wired into the UI** — the headline number, on held-out failures.

**Phase 2 — lower-confidence outcome breadth for the metabolism/phenotypic gap (IMPLEMENTED)**
- [x] Reactive-metabolite structural alerts plus hepatotox/mito read-across in `app/outcome_modules.py`.
- [x] Fuse the off-target core + outcome modules into one UI/report, each item tagged by tier + confidence; liver/mito shown as *model-predicted*, not mechanism-linked, and never merged into the validated off-target headline.

**Phase 3 — Rigor & hardening (post-MVP)**
- [ ] Per-target AUC beyond the 4 measured; swap class-Tanimoto → trained target-prediction (SEA / per-target QSAR); metabolite site-of-metabolism prediction.
- [ ] Calibration; matched No-concern comparator; scale curation to ~200 drugs (~2 person-weeks).

**Phase 4 — Stretch**
- [ ] Learned-embedding method (only with demonstrated scaffold-split lift); off-target docking panel; AACT `why_stopped` mining; 3D O3A/ROCS test.

---

## 6. Demo / pitch (investor-facing)

**Claim:** *"For a candidate that's going to die on a safety liability, we run the assay that kills it **first** — cutting the experiments, time, and cost you burn before the go/no-go decision. We do it by linking the molecule to clinical failures by *shared mechanism*, including liabilities no structure search would find."* The headline is a **number: on real historical failures we cut assays until culprit from ~11 to ~4.**

Framed honestly as a **retrospective leave-one-out demo** — hide a known failure's label, query it, show the system recovers the liability *and* reorders the assay plan. Two levels:

**Primary demo — a buried liability (the real value): rimonabant (or pergolide).**
1. Paste **rimonabant** as the "candidate," its tox label hidden (leave-one-out; its CB1 partner removed too).
2. *Screening against N toxicity-causal clinical failures vs. matched No-concern drugs…*
3. **The value moment:** *"A standard secondary-pharmacology panel front-loads the cardiac channels; the CB1 counter-screen that would have caught rimonabant's psychiatric liability sits at rank 15. Our tool flags CB1 at **z = +10.8** and moves it to **rank 1** — a *no-go*, because rimonabant resembles drugs withdrawn for psychiatric toxicity. That's ~11 fewer assays before the kill."* (For pergolide the same story runs on 5-HT2B valvulopathy, default rank 11 → rank 1 — historically resonant, the fen-phen mechanism.)
4. **Assays-until-culprit gain, shown as the headline:** default plan rank → our rank, per drug, with the aggregate **11.3 → 3.8**.
5. Each item carries the driving failed drug + tier + citation.

**Validation anchor (not the headline) — terfenadine.** Show the clean recovery (z = +6.4 on hERG even with terfenadine + cisapride removed) as proof the engine works — but state plainly that hERG is *already* front-loaded, so this is a validation of the method, not a reordering win. **Do not hang the value claim on hERG.**

**Trust signals (raise them before a judge does):**
6. Paste an out-of-domain molecule (cyclosporine) → the tool **abstains** ("outside small-molecule applicability domain").
7. Paste a candidate that closely resembles a known failure → the tool flags **known analog** ("cheap similarity already catches this — our engine adds most for novel chemotypes"). This *pre-empts* the ablation critique instead of hiding from it.
8. Mechanism-edge network; click a node to inspect a failed drug.

**Value prop:** broad safety panels are expensive and slow — *"move the assay that will kill the program to the front"* saves cost, time, animals, and reaches go/no-go earlier. Credible because it's grounded in real withdrawals, measured on 20 historical failures against a defensible baseline, ablated against dumb similarity, and honest about where it adds nothing (hERG) and where it's blind (metabolite tox).

---

## 7. Risks & honest framing (raise these before a judge does)

1. **Curated failure-label provenance.** Withdrawn ≠ toxicity-caused. → mandatory `failure_reason` + `safety_related` curation before treating a drug as a toxicity-causal failure (§4 schema).
2. **Small curated failure class** (~90–110 in ClinTox). → lean on organ-liability labels (DILIrank/DICTrank ~1.3k) + off-target class matching for scoring; curated failures provide the clinical-grounding badge.
3. **"Safe" comparator is a misnomer.** → matched **No-concern** class, not "any approved" (§3e).
4. **Score semantics.** Present enrichment (`N×`) or 0–100 priority, **not** probability-looking decimals unless calibrated.
5. **No value on hERG.** Measured: default panels already front-load it, and we're sometimes worse. → frame the value as *buried* liabilities (5-HT2B, CB1, 5-HT3); use hERG only as a validation anchor.
6. **The engine only beats dumb similarity on novel chemotypes.** Measured ablation: naive 2D ties/wins when an analog exists; off-target class matching wins (naive 0/20 → class matching 12/20) only under strict hold-out. → state the claim conditionally; ship the **known analog flag** so the tool is honest about when it adds value; measured off-targets for reference drugs, predicted only for the candidate.
7. **Off-target prediction is itself a model** — worst for novel scaffolds (AUC 0.667 diverse-isolated holdout vs 0.913 pooled). → report confidence + applicability-domain flag.
8. **Some clinical killers are NOT off-target binding** — reactive/CYP metabolites, idiosyncratic mito tox, dose/exposure. Measured: our only misses (fenfluramine, dexfenfluramine, sibutramine) are metabolite-active. → these are handled by **lower-confidence outcome modules**; liver/mito is a *designed module boundary*, but the outcome modules are less validated than the off-target core — don't oversell them.
9. **Generalization.** Never quote a random-split number; report scaffold-split CV + AD abstention.
10. **Narrow addressable population.** Measured in a blind hold-out (§4): drugs withdrawn *specifically for a reachable off-target* are a minority — the objective withdrawn population is dominated by off-panel hepato/hematologic/carcinogenic tox and on-target pharmacology. → the value metric (assays until culprit) applies to *off-target-mediated* failures, not "all safety failures"; size the market/claim on that slice, not on total withdrawals.
11. **Proprietary gap.** Cleanest discontinued-for-safety data (Pharmaprojects/Citeline, Pharmapendium) is paywalled — note as gold standard if licensed.
12. **Licenses:** DrugBank / SIDER / OFFSIDES non-commercial — flag for any commercialization.

---

## 8. Tools
**Core:** RDKit (ECFP/FCFP, Gobbi pharmacophore, USRCAT 3D, conformers via ETKDGv3, descriptors, substructure), **FPSim2** (2D retrieval), scikit-learn (scaffold-split CV, calibration, stacking), an LLM for the agent report. **Env managed with `uv`** (venv at `.venv`, rdkit 2026.03.3; `uv pip install --python .venv …`).
**Off-target class matching:** **Bowes 2012 ~44-target panel** (aligned to SAFETYscan47) as the vocabulary; **ChEMBL** actives per target (pChEMBL ≥ 6, IC50/Ki/Kd) via REST; SEA / ChEMBL target-prediction / PIDGIN and per-target QSAR (hERG has strong public models) for the production upgrade; optional docking.
**Baselines / cross-checks:** naive 2D nearest-failed-drug (ablation), plain ECFP, ProTox 3.0 (driver in `experiments/derisk/protox/`), ADMET-AI, EPA GenRA.
**Structure resolution:** PubChem PUG REST (`SMILES`/`ConnectivitySMILES` properties), ChEMBL.
**Only-if-lift-shown:** MolFormer/ChemBERTa/Uni-Mol as an added method, gated by a scaffold-split holdout.

---

## 9. Key sources
- ClinTox — MoleculeNet arXiv 1703.00564; TDC `single_pred/Tox`
- WITHDRAWN 2.0 — NAR 2024;52:D1503 (v1: NAR 2016;44:D1080)
- Probes & Drugs — probes-drugs.org/compoundsets · AACT — aact.ctti-clinicaltrials.org
- DILIrank / DILIst — FDA LTKB · DICTrank — Drug Discov Today DOI 10.1016/j.drudis.2023.103770
- Tox21 / ToxCast — tox21.gov; EPA figshare invitrodb · SIDER / OFFSIDES / TWOSIDES — nsides.io; Tatonetti Lab
- **Closest prior art** — Albert, Skinner et al., *ACS Pharmacol Transl Sci* 2025 — **DOI 10.1021/acsptsci.5c00452, PMID 40969870** (SAFETYscan47; 52 safety-discontinued compounds)
- **Prior art to differentiate** — AbbVie OTSA; ToxEvaluator (Pfizer 2016); Deaton 2019; Zatorski 2025 / Xue 2013 (see [`novelty-search-log.md`](novelty-search-log.md))
- Secondary-pharmacology safety panel vocabulary — Bowes et al., *Nat Rev Drug Discov* 2012;11:909
- Predictors — ADMET-AI; ProTox 3.0 (NAR 2024 W513); ADMETlab 3.0 (NAR 2024 W422) · Read-across — EPA GenRA

**Verified demo mechanisms — citations:**
- CB1 / rimonabant psychiatric — Christensen *Lancet* 2007 (PMID 18029027); taranabant — Aronne 2010
- 5-HT2B valvulopathy — Setola et al., *Mol Pharmacol* 2003 (PMID 12761331); pergolide — Zanettini *NEJM* 2007 (PMID 17202454); fenfluramine — Connolly *NEJM* 1997 (PMID 9271479); benfluorex — Tribouilloy 2011 (PMID 21981882)
- 5-HT3 / alosetron ischemic colitis — FDA 2000 withdrawal/restriction
- hERG shared determinants (terfenadine + cisapride) — Saxena et al., PMC2845965
- Dual BSEP+mito DILI (troglitazone + nefazodone) — Aleo et al., *Hepatology* 2014;60:1015
- SMILES/IDs from PubChem & ChEMBL; ECFP4/FCFP Tanimoto + z-scores computed locally (RDKit, Morgan r=2, 2048-bit; scripts in `experiments/`)

---

## Appendix — de-risking experiments (reproducible)
All in [`experiments/`](experiments/) (env: `uv` venv, RDKit 2026.03.3). Each `experiments/derisk/*/` has a `FINDINGS.md`.

| Experiment | Folder | Key result |
|---|---|---|
| 2D / 3D / off-target class-matching feasibility | `experiments/` | 2D fails (0.05–0.19); class matching recovers (z=+6.4) |
| Specificity + cross-target + scaffold + coverage | `experiments/` | AUC 0.894/0.95/0.92/0.88; scaffold 0.913/0.667; 37/39 |
| **Assay-recovery n=20 (decisive)** | `derisk/expanded/` | buried assays until culprit **11.3→3.8**; 6 non-obvious wins |
| **Off-target class matching vs 2D ablation** | `derisk/ablation/` | class matching necessary only for novel chemotypes (naive 0/20 strict leave-one-out) |
| **ProTox competitor** | `derisk/protox/` | partial cardiac; can't localize buried; not redundant |
| Abstain / applicability domain | `derisk/abstain/` | descriptor-box rule 0/31, 0/11; NN rule fails |
| Metabolite handling | `derisk/metabolite/` | MAX(parent,metab); fenfluramine +0.35→+1.22 |
| **Blind hold-out (external selection + labels)** | `derisk/blind_holdout/` | magnitude replicates (4.7 vs 11.4); addressable class narrow — off-panel/on-target dominate |
| Novelty search log | `novelty-search-log.md` | triple-combination white space survives |
