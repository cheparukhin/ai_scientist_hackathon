# Toxicity Risk → Assay Recommender

**A retrieval tool that scores a small-molecule candidate's toxicity risk by linking it to drugs that failed in the clinic (or were withdrawn) for toxicity, and recommends which in-vitro assays to run first.**

> Positioning in one line: *"We connect your candidate to clinically-failed drugs it shares a toxicity **mechanism** with — including ones no 2D-structure search would find — and rank which in-vitro assays would catch the implicated liability first. It prioritizes experiments; it does not declare a molecule safe."*

---

## 0. Project context & status (read first)

- **Event:** a bio × AI hackathon. Deliverable is an **MVP + investor/judge-facing demo** (see §6), not a validated product.
- **How we got here:** chosen over an alternative — a *cross-specificity / immunogenicity risk tool for biologics* (antibody candidates, e.g. cetuximab/α-gal) — because the small-molecule tox tool is **more tractable**: core datasets are downloadable, and the MVP is retrieval + light reasoning. The biologics idea is a fallback.
- **Key pivot (from our own feasibility test — see §4):** we tested whether cheap **2D structural similarity** actually connects mechanistically-related tox failures. **It does not** (7 verified withdrawn-drug pairs, ECFP4 Tanimoto 0.05–0.19). The shared liability is an off-target *binding* phenomenon 2D structure can't see. So the retrieval engine is **mechanism / off-target similarity**, with 2D fingerprints demoted to a fast "obvious-analog" pre-filter.
- **Repo state:** the repo contains this doc + an **[`experiments/`](experiments/)** folder with reproducible feasibility scripts (2D / 3D / R4) + results (RDKit 2026.03.3, seeds, SMILES). No product code yet — the pipeline (§4) is greenfield.
- **Guiding principle:** a **portfolio, tiered by confidence** (§4). The validated off-target engine (M1/R4) is the demo headline and covers *broad* off-target liabilities (cardiac/CNS/GI/endocrine/immune — 37/39 panel targets); wrapped outcome-models (M2) cover the liver/mito gap; every recommendation is tagged **mechanism-linked / model-predicted / abstain**. The 2D floor is a safety-net, not the pitch.
- **Parameters still TBD:** timeframe/deadline; team size & strengths; explicit judging criteria; target LLM + infra.
- **Data caveat:** dataset sizes are **last-published figures** and several conflict across mirrors (e.g. ClinTox is cited as **1,484** on TDC and **1,491** in the MoleculeNet paper). **Re-confirm every count and DOI against the primary source before quoting in a pitch.**

---

## 1. The idea

- **Input:** a small-molecule candidate (SMILES; optionally a 3D conformer).
- **Reference database:** drugs annotated with **why they failed** — prioritizing those that **failed clinical trials or were withdrawn specifically for a toxicity**, each tagged with the **organ, phenotype, and mechanism** and a provenance link.
- **Output:**
  1. A **per-modality risk profile** expressed as **enrichment vs. a matched lower-liability comparator** — e.g. *"liver-liability enrichment 4.1× (95% CI …); cardiac 1.2×"* — or a **0–100 priority score**, **not** a probability of harm (see §4 Layer 3).
  2. A **ranked in-vitro assay panel** — e.g. *"prioritize hERG + iPSC-cardiomyocyte (candidate shares predicted hERG liability with terfenadine/cisapride, both withdrawn for QT); then hepatocyte spheroid."*
  3. An **evidence trail** — the linked failed drugs, the organ/mechanism they failed by, whether the link is structural or mechanistic, and provenance.

It is a **decision-support / experiment-prioritization** tool, not a binary safety classifier — more honest (structure/mechanism proxies miss metabolism, dose, exposure) and easier to validate and adopt.

---

## 2. Why it's novel (white space + prior art to differentiate from)

The two halves each exist separately; the **combination is unoccupied**:

- **Endpoint predictors** — ProTox 3.0, admetSAR 3.0 / ADMETlab 3.0, ADMET-AI, DeepTox, DEREK. Predict *tox probabilities / structural alerts*, **not "which assay to run,"** and are **not grounded in clinical-trial failures**.
- **Read-across engines** — EPA **GenRA**, OECD QSAR Toolbox, AMBIT. Analog-based tox inference (closest cousin), but reference sets are **environmental chemicals**, and output isn't a prioritized clinical assay panel.
- **Closest paper** — Albert, Skinner et al., *"Comprehensive Analysis of Clinically Discontinued Compounds Using an In Vitro Secondary Pharmacology Panel…"* (**ACS Pharmacol Transl Sci 2025, DOI 10.1021/acsptsci.5c00452**): **52 compounds discontinued 2001–2021 for safety** run through the **SAFETYscan47** (47-target) panel. It screens **everything** and argues the primary target doesn't predict off-target risk — it does **not** use similarity to *prioritize a subset* of assays. **That is our niche.**

> **Novelty claim — needs a documented search log, don't overstate.** A *preliminary* search found no public tool/paper combining (clinical-failure grounding) × (similarity/mechanism retrieval) × (specific assay recommendation). Before the pitch, record the exact **databases, dates, query strings, and inclusion criteria** so "we found nothing equivalent" is reproducible rather than an assertion.

**Our two defensible novelty claims:**
1. **Reference set grounded in *toxicity-causal* clinical failures** (not assay-positives, not mere market withdrawals).
2. **Output = a specific, ranked assay panel**, made credible by organ+mechanism-resolved labels and a validation protocol.

> **Be precise about T1 vs T2 in the pitch:** the *clinical-failure grounding* (T1) is what's novel, but the *scoring backbone* is largely **T2 organ-liability labels** (DILIrank/DICTrank — marketed drugs). Don't imply the whole score comes from toxicity-causal clinical failures — T1 is the grounding badge, T2 does the statistical lifting.

Genuine white space = the **linkage layer**: candidate → the mechanism/off-target it shares with clinically-failed drugs → the assay that would have caught them.

---

## 3. Datasets

### 3a. Toxicity-causal failure set (the novelty; assemble + curate, don't just merge)
No public dataset cleanly says "these N molecules died in the clinic *because of toxicity*, here are their SMILES + organ + mechanism." Build it — and **curate the *reason*, don't assume it**:

| Dataset | Content | Size | Access | Role |
|---|---|---|---|---|
| **ClinTox** (MoleculeNet / TDC) | FDA-approved vs. **failed clinical trials for toxicity** | **1,484 (TDC) / 1,491 (paper)**; tox-failure positives only ~90–110 | Free — `deepchem`, TDC, HF `zpn/clintox` | Start here. Positives mined from AACT — our exact idea, v0. Class-imbalanced. |
| **WITHDRAWN 2.0** (Charité) | Withdrawn/discontinued drugs + withdrawal reasons | ~578 (v1) → ~631–635 (2.0) | Free — cheminfo.charite.de/withdrawn; NAR 2024 52:D1503 | ⚠️ **Includes *non-safety* withdrawals** (commercial, efficacy, superseded). Must filter on `safety_related` + organ before calling anything T1. |
| **Probes & Drugs** — withdrawn set | Withdrawn drugs pre-joined to ChEMBL/DrugBank | ~633 (verify release/date) | Free — probes-drugs.org/compoundsets | De-dupe / canonical-SMILES join. |
| **DrugBank** (withdrawn flag) | `withdrawn` category, structures | ~250–300 (verify) | Free academic; **non-commercial** | Structure + flag; weak on reason. |
| **AACT** (ClinicalTrials.gov dump) | `why_stopped` free-text, `TERMINATED/…`, interventions | full CT.gov | Free — aact.ctti-clinicaltrials.org | **Stretch.** NLP-classify tox terminations → resolve names to SMILES (~30–60%). Demo a slice; don't rebuild wholesale. |

### 3b. Organ + mechanism-resolved toxicity labels (power the "which assay" step)
FDA-curated *ranked* lists — highest-quality organ-resolved labels; backbone of scoring **and** of the matched comparator (§3e).

| Dataset | Content | Size | Assay it points to |
|---|---|---|---|
| **DILIrank** (FDA LTKB) | **liver injury**, 4 classes incl. an explicit **No-DILI-concern** class | ~1,036 (v1); **2.0 = 1,336** | Hepatocyte/spheroid, BSEP inhibition, mito-tox |
| **DILIst** (FDA) | binary hepatotoxic (broader) | ~1,279 | same |
| **DICTrank** (FDA, 2023) | **cardiotoxicity**, 4 concern levels incl. **No-DICT-concern** | ~1,318 | hERG, iPSC-cardiomyocyte, Nav1.5/Cav1.2 · DOI 10.1016/j.drudis.2023.103770 |
| **hERG datasets** | hERG channel block | TDC 648; DeepHIT ~14k; UnihERG ~20k | hERG patch-clamp / electrophysiology |
| **SIDER** | side effects from labels | 1,556 drugs | phenotype → organ → assay |
| **OFFSIDES / TWOSIDES** | FAERS-derived side effects / DDIs | 438,801 / 868k assocs | pharmacovigilance layer |
| **LiverTox** (NIH) | narrative hepatotox monographs | ~1,200 agents | qualitative liver context (text) |

> The **No-concern** classes in DILIrank/DICTrank are the best public **lower-liability comparators** (see §3e) — much better than "any approved drug."

### 3c. Assay / mechanism datasets (map mechanism → a named assay)
| Dataset | Content | Size | Role |
|---|---|---|---|
| **Tox21** | 12 pathway assays | ~7,800 × 12 | each pathway *is* an assay |
| **ToxCast** (EPA) | HTS in-vitro endpoints | ~800–1,400 endpoints × >4,200 chem (verify per invitrodb release) | richest compound → assay-endpoint matrix |
| **ChEMBL** | measured bioactivities incl. safety targets (hERG, 5-HT2B, BSEP…) | ~2.4M cpds | **the mechanism layer** — reference drugs' *known* off-targets; MCP available |
| TG-GATEs / DrugMatrix / LINCS | toxicogenomics / histopath / cell response | ~170 / ~600 / >20k | mechanism enrichment (stretch) |

### 3d. Structure sources
**PubChem** (name→structure resolution), **ChEMBL** (CC-BY-SA; structures + bioactivity; MCP), **DrugBank** (annotations; non-commercial).

### 3e. Matched lower-liability comparator (the denominator — do not call it "safe")
The class-conditional score needs a denominator. **"Approved" ≠ "safe"** — approved drugs carry organ liabilities too. Use **per-organ No-concern classes** (DILIrank *No-DILI-concern*, DICTrank *No-DICT-concern*) as the comparator, **matched to the toxic set by scaffold / physchem / route** where possible, so the enrichment reflects *toxicity signal*, not just drug-likeness. Without a matched comparator, everything looks risky.

---

## 4. Architecture

### Empirical reality check (we ran this — it drove the design)
We assembled 7 pairs of drugs that were **both withdrawn/failed for the same organ toxicity via the same documented mechanism**, then computed similarity. **2D fingerprints miss all of them:**

| Pair | Shared mechanism → organ | ECFP4 | FCFP4 | Clinical status (precise — verify per row) |
|---|---|---|---|---|
| **Terfenadine × Cisapride** | hERG block (Tyr652/Phe656) → QT/TdP | 0.19 | 0.15 | both withdrawn (US 1998 / 2000) |
| Astemizole × Thioridazine | hERG → QT/TdP | 0.17 | 0.25 | astemizole withdrawn 1999; thioridazine brand withdrawn 2005 (QT boxed 2000) |
| **Pergolide × Fenfluramine** | 5-HT2B agonism → valve fibrosis | 0.11 | 0.14 | both withdrawn (US 2007 / 1997) |
| Cabergoline × Benfluorex | 5-HT2B (norfenfluramine) → valvulopathy | 0.11 | 0.15 | cabergoline **restricted/boxed, not withdrawn** (still marketed low-dose); benfluorex withdrawn EU 2009 (Mediator) |
| **Troglitazone × Nefazodone** | BSEP + mitochondrial → DILI | 0.13 | — | troglitazone withdrawn 2000; nefazodone brand discontinued ~2004 (hepatotox boxed) — both DILIrank-Most |
| Troglitazone × Sitaxentan | BSEP inhibition → cholestatic DILI | 0.11 | — | both withdrawn (2000 / 2010) |
| Fialuridine × Perhexiline | mitochondrial → hepatotox | 0.05 | — | FIAU **never approved** (Phase 2 halted 1993, 5 deaths); perhexiline withdrawn most markets, TDM-only AU/NZ |

Feature/pharmacophore 2D fingerprints (FCFP, Gobbi) gave only **1.3–1.5× lift** (absolute sim still ~0.15–0.25); sanity checks passed (true analogs stayed high: pergolide↔cabergoline 0.41). **Conclusion:** the signal linking these is **shared binding target**, which needs a 3D and/or off-target layer — *not* 2D structure. This is the single most important design fact.

**3D shape test (USRCAT) — done, mixed result.** 30 conformers/molecule (ETKDGv3 + MMFF), each pair's true partner ranked against the 16-compound library. USRCAT improved partner rank in **4/7** pairs and lifted top-3 hits only 3/7 → 4/7. One clean rescue: **fialuridine–perhexiline** (mito) — ECFP4 buried it at rank 11 (Tanimoto 0.05) → USRCAT **rank 1**. But it *worsened* 3 pairs (astemizole–thioridazine 2→9, troglitazone–nefazodone 2→9) because **global molecular shape/size dominates the descriptor, not the target-relevant toxicophore**; absolute sibling scores (0.12–0.23) stay non-discriminative. **Takeaway: 3D shape is a weak *supplement* (helps shape-driven cases), not a dependable rescue → R4 (off-target/mechanism) remains the engine.** A stricter 3D test — O3A/ROCS alignment or explicit pharmacophore-fitting, which target the 3-point toxicophore geometry USRCAT ignores — is the remaining untested 3D option.

**R4 mechanism-linkage test (the fix) — done, decisive.** We prototyped the off-target layer as a lightweight proxy: for each safety target, fetch its known actives from ChEMBL (hERG 1,483; 5-HT2B 1,214; Nav1.5 258 — **BSEP yielded only 5, so it fell back**), and score a query by mean top-5 Tanimoto to that *class* (leave-one-out: query + named partners removed first). Result: **terfenadine → hERG = 0.57 (z = +6.4 vs background)** vs pairwise terfenadine–cisapride = 0.19 — and cisapride wasn't even in the fetched hERG set, so the recovery is driven entirely by *other* hERG ligands. Across the 6 query drugs (each scored against hERG, 5-HT2B, Nav1.5): the correct mechanism was the **single highest-scoring target for 4 of 6**, and scored **above the non-binder background for all 6** — i.e. always detectable, and the outright top pick two-thirds of the time. The 2 non-top cases are honest and explained: fenfluramine (its active is the metabolite norfenfluramine, so the parent scores low) and thioridazine (genuinely both a hERG blocker *and* serotonergic, so it ties rather than errs). **Conclusion: class-aggregation recovers the mechanism pairwise buries — R4 is validated and cheap.** Full scripts + data in [`experiments/`](experiments/).

### De-risking results — four cheap checks before building (all run; scripts in [`experiments/`](experiments/))
Before committing to a build we stress-tested the load-bearing assumptions:

| Check | Question | Result | Verdict |
|---|---|---|---|
| **Specificity** (hERG) | Does the score separate blockers from non-blockers, or flag everything? | **ROC-AUC 0.894** — blockers (median 0.66) vs. assay-confirmed non-blockers (0.28); 19/19 common drugs scored low | ✅ discriminates; no false-positive flood |
| **Scaffold generalization** | Does it survive *novel* chemotypes (the real use case)? | Pooled leave-one-scaffold-out AUC **0.913**; a diverse single-ring holdout collapses to **0.667** | ⚠️ real, but degrades for truly novel, isolated chemotypes |
| **Panel coverage** | How much of the safety panel can R4 serve? | **37/39 targets serviceable (34 rich)** — all aminergic GPCRs, transporters, enzymes, nuclear receptors. Gap: **BSEP 5 actives; mitochondrial = no single target** | ✅ broad off-target reach; ❌ liver is a data desert |
| **Phase-0 curation** | Cheap merge or manual sink? | SMILES 20/20 auto; correct mechanism auto-derivable **7/20** (all the off-target cases); `failure_reason`/`tier` ~0% structured | ⚠️ hybrid: cheap for off-target drugs, ~2 person-weeks curation for the rest |

**What this means for scope.** R4's off-target engine is **broad, not narrow** — 37/39 secondary-pharmacology targets spanning cardiac, CNS, GI, endocrine, and immune liabilities, and it genuinely discriminates (AUC 0.89 on hERG, the one target measured end-to-end). The honest boundary is **off-target-mediated (R4's domain) vs. metabolism/phenotypic**: liver DILI (BSEP/mito) and reactive-metabolite/idiosyncratic tox are *not* off-target-binding events, so no similarity method reaches them — they need outcome-trained models (the portfolio below). *Caveat:* coverage ≠ measured accuracy — only hERG has a validated AUC; per-target discrimination across the other 36 rich targets is expected (aminergic GPCRs are QSAR-tractable) but unmeasured.

### The tool is a *portfolio*, not one engine
Breadth comes from the right method per domain, fused into one prioritized assay panel, each item tagged with an **evidence tier**:

| Module | Domains | Method | Evidence tier (shown in UI) |
|---|---|---|---|
| **M1 — off-target linkage (R4, the validated core)** | cardiac, CNS, GI, endocrine, immune (37 targets) | similarity to known actives → linked clinical failures | **mechanism-linked** (strongest) |
| **M2 — outcome models** | liver DILI, mitochondrial, reactive-metabolite | QSAR on DILIrank, Tox21-mito model, structural alerts; wrap ADMET-AI / ProTox | **model-predicted / alert-based** |
| **Fusion** | all | rank + tag every recommendation by evidence tier & confidence; **abstain** outside domain | — |

The novelty is the **fusion + clinical-failure grounding + assay triage + honest evidence-tiering**, not any single predictor. Coverage flags are first-class output: *covered mechanism / weak coverage / abstain*.

### Pipeline
```
             ┌─────────────────────────────────────────────┐
  SMILES  →  │ 1. Multi-rung similarity / linkage           │
             │   R1  2D fingerprint (ECFP, FPSim2)  [floor] │  ← fast pre-filter, obvious analogs
             │   R2  feature/pharmacophore fp (FCFP)        │
             │   R3  3D shape + pharmacophore (USRCAT)      │  ← test-gated
             │   R4  off-target / mechanism similarity      │  ← THE differentiator
             │       (predicted 2ary-pharmacology profile)  │
             └───────────────────┬─────────────────────────┘
                                 ▼
             ┌─────────────────────────────────────────────┐
             │ 2. Link to failed-drug reference DB          │
             │   • candidate ⇄ drugs w/ shared mechanism    │
             │   • carry organ + phenotype + tier + source  │
             └───────────────────┬─────────────────────────┘
                                 ▼
             ┌─────────────────────────────────────────────┐
             │ 3. Per-modality scoring (class-conditional)  │
             │   • enrichment vs matched No-concern class   │
             │   • report as N× enrichment or 0–100 priority │
             │   • applicability-domain flag                │
             └───────────────────┬─────────────────────────┘
                                 ▼
             ┌─────────────────────────────────────────────┐
             │ 4. organ/mechanism → assay mapping           │
             │   → ranked assay panel w/ driving evidence   │
             └───────────────────┬─────────────────────────┘
                                 ▼
             ┌─────────────────────────────────────────────┐
             │ 5. Explainability + Agent + Visualization    │
             └─────────────────────────────────────────────┘
```

### Reference DB schema (split labels — avoid circular recommendations)
Keep organ, phenotype, mechanism, and assay as **separate fields** so the recommended assay is *derived*, not baked into, the label:
```
compound_id | canonical_SMILES | source(s) | provenance/evidence_link
failure_reason (free text) | safety_related (bool) | fail_stage | label_tier
organ_system | tox_phenotype | mechanism (e.g. hERG block, BSEP inhib, 5-HT2B agonism)
known_offtargets (from ChEMBL) | assay_endpoint (observed) | recommended_assay (derived)
```
- **`safety_related` + `failure_reason` are mandatory** — WITHDRAWN/discontinued includes non-safety reasons; never tier a compound T1 without a curated toxicity-causal reason + provenance.
- **Label tiers** (never blend silently): **T1** curated toxicity-causal failure (safety_related=true, organ+mechanism, evidence link); **T2** organ-severity-ranked marketed (DILIrank/DICTrank Most); **T3** noisier flags (SIDER-derived, broad ClinTox).

### Layer 1 — Multi-rung similarity / linkage
Each rung is interpretable and maps to a known driver; the **mechanism rung (R4) is the engine**, per the reality check.

| Rung | Captures | Tool | Role |
|---|---|---|---|
| **R1 2D fingerprint** (ECFP4 Tanimoto) | shared 2D substructure | RDKit + **FPSim2** | fast floor; catches *obvious* analogs only |
| **R2 feature/pharmacophore fp** (FCFP, ErG/Gobbi) | pharmacophoric features (topological) | RDKit | marginal lift; cheap |
| **R3 3D shape** (USRCAT; O3A/ROCS = stricter, untested) | conformer shape (global, not toxicophore geometry) | RDKit USRCAT | **tested → weak supplement**: rescued a shape-driven mito pair (rank 11→1) but hurt others; not dependable |
| **R4 off-target / mechanism** | shared (predicted) protein liabilities: hERG, 5-HT2B, BSEP, mito… | target-prediction (SEA, ChEMBL target-prediction, PIDGIN) or a small **docking panel**; matched to reference drugs' ChEMBL off-targets | **the differentiator** — this is what links the verified pairs |

**Layer 1 — R4 detail (the safety-panel pipeline).** The mechanism vocabulary is a **defined secondary-pharmacology safety panel**, not an open-ended target list: the **Bowes et al. ~44-target minimal panel** (*Nat Rev Drug Discov* 2012;11:909), aligned to the **SAFETYscan47** panel in the discontinued-compounds paper (*ACS Pharmacol Transl Sci* 2025, DOI 10.1021/acsptsci.5c00452). Each off-target (hERG, 5-HT2B, BSEP, Nav1.5, mitochondrial complexes, …) maps to an **organ** and a **named assay** (see Layer 4). Pipeline:
- **(a) Reference side** — annotate each clinically-failed drug with its **measured** panel activity from ChEMBL (hERG/5-HT2B/BSEP IC50 etc.) + organ/mechanism/provenance. These are facts, not predictions.
- **(b) Candidate side** — score the candidate against the **same panel** via off-target/target-prediction (SEA, ChEMBL target-prediction, PIDGIN) and/or per-target QSAR (hERG has excellent public models — see §3b).
- **(c) Link** — candidate's predicted panel-hits → failed drugs sharing them → organ → assay.

**Concrete MVP path (validated — commit to this one).** For each safety-panel target, fetch known actives from **ChEMBL** (pChEMBL ≥ 6; IC50/Ki/Kd); score the candidate by **mean top-5 Tanimoto to the class**, reported as a **z-score vs. a background set of non-binder drugs** (prototype: terfenadine→hERG z = +6.4). **Fallbacks:** if a target has too few actives (BSEP returned only 5 → unusable), fall back to that mechanism's direct reference-drug annotations or the 2D floor, and flag low target-coverage. **Production upgrade:** swap raw class-Tanimoto for a trained target-prediction model (SEA / per-target QSAR) and handle metabolite-active drugs explicitly (the fenfluramine miss). Start the panel with **hERG and 5-HT2B / aminergic GPCRs** — best-covered, and the verified-pair mechanisms.

**Why R4 works where pairwise similarity failed.** Terfenadine and cisapride aren't similar to *each other* (ECFP4 0.19) — but each is similar to *some* member of the large, structurally diverse set of known hERG ligands. Predicting "binds hERG" by comparing to the **whole target class** is far more sensitive than comparing to one sibling. R4 wins by **aggregating over the target's entire diverse ligand set** — exactly the signal a single pairwise comparison throws away. *The off-target predictors are the means, not the product:* the novelty is the clinical-failure grounding + linkage + assay recommendation, not the predictors themselves.

**On learned embeddings:** off-the-shelf encoders (MolFormer/ChemBERTa/Uni-Mol) embed *general chemistry, not toxicity*; unproven lift, unexplainable. Include only if they beat R1–R4 on a scaffold-split holdout. Until then, cut them.

### Layer 3 — Per-modality scoring (class-conditional; enrichment, not a probability)
Similarity to a toxic drug means nothing unless the candidate is *more* linked to toxic drugs than to the **matched No-concern comparator** (§3e). Per organ *m*:

```
enrichment_m  =  linkage(candidate, T1/T2 failures in organ m)
                 ─────────────────────────────────────────────
                 linkage(candidate, matched No-concern comparator)
```

Report as **"4.1× liver-liability enrichment"** or a **0–100 priority score** — never `0.71`-style decimals that read as probability of harm, unless calibrated. Combine rungs by transparent rank-aggregation/logistic stack (not hidden weights). Attach: driving links (drug, organ, mechanism, tier, provenance) + **applicability-domain flag**.

### Layer 4 — organ/mechanism → assay mapping (derived, hand-authored)
| Mechanism / organ | Recommended assays (priority order) |
|---|---|
| hERG block / cardiac (DICTrank) | hERG patch-clamp → iPSC-cardiomyocyte → Nav1.5/Cav1.2 |
| BSEP inhibition / cholestatic liver | BSEP (bile-salt-export-pump) inhibition → hepatocyte spheroid |
| Mitochondrial / hepatocellular liver | mito assays (Glu/Gal, Seahorse) → hepatocyte spheroid |
| 5-HT2B agonism / valvulopathy | 5-HT2B binding/functional counter-screen (per Setola 2003) |
| Reactive metabolite | GSH-trapping / covalent-binding; CYP-dependent cytotox |
| Tox21/ToxCast pathway hit | the specific named HTS assay for that pathway |

### Layer 5 — Explainability, agent, visualization
- **Attribution:** name the *shared mechanism/off-target* and (where structural) the driving substructure — fingerprint-bit→atom mapping for R1/R2, docked pose or shared pharmacophore for R3/R4. Defends against "isn't this just fingerprint similarity?"
- **LLM agent** turns retrieved evidence into a report (retrieval + reasoning; no training needed for v1).
- **Visualization** — candidate at center, linked failed drugs around it, **edges labeled by mechanism** and colored by organ; click a node → why it failed, mechanism, literature, suggested assay.

### Rigor: validation, calibration & applicability domain
- **Scaffold-split CV (never random)** — random splits leak analogs and inflate everything; Bemis–Murcko split tests generalization to novel chemotypes.
- **Per-organ discrimination:** ROC-AUC and **enrichment factor / precision@k** (triage-relevant).
- **Assay-recovery (tests the actual product promise):** for held-out failures with a known mechanism, does the **correct assay appear in the top-1 / top-3** of the recommended panel? This validates the *"which assay first"* claim directly — organ discrimination alone doesn't.
- **Calibration:** does 4× enrichment map to a real elevated failure rate? Report a calibration curve; ranked-but-uncalibrated is fine for *prioritization* if not sold as P(harm).
- **Baselines to beat:** plain ECFP (does R3/R4 add signal?) and ≥1 published predictor (ProTox / ADMET-AI) on the same holdout — demonstrate **lift**.
- **Applicability domain (OECD):** flag and **abstain** outside the reference set's chemical/target space.
- **Ablation:** report each rung's contribution — "R4 is what caught this," measured not asserted.
- **Honest bounds:** structure/mechanism proxies miss metabolism, reactive metabolites, dose/exposure; T1 positive class is small (lean on T2 for scoring, T1 as the clinical-grounding badge).

---

## 5. Roadmap (revised after de-risking)

The de-risking (Phase 0) is **done** and reshaped the plan: build a **tiered portfolio** — ship the validated off-target core (M1) first, add the metabolism/phenotypic modules (M2) for breadth, then harden.

**Phase 0 — De-risking (DONE ✅ — see §4 + [`experiments/`](experiments/))**
- [x] 2D similarity fails on mechanistic pairs (ECFP4 0.05–0.19); 3D/USRCAT a weak supplement (helped 4/7, hurt 3/7).
- [x] R4 off-target linkage recovers the mechanism (terfenadine→hERG z = +6.4; top target 4/6, above background 6/6).
- [x] **Specificity AUC 0.894** (blockers vs non-blockers); scaffold-generalization 0.913 pooled, but 0.667 for truly novel isolated chemotypes.
- [x] **Coverage 37/39 panel targets serviceable**; liver (BSEP/mito) is the gap.
- [x] Phase-0 curation is a hybrid (cheap for off-target drugs; ~2 person-weeks for mechanism/reason/tier).
- → **Verdict: viable, broad within off-target space, honest liver/mito gap → build as a tiered portfolio.**

**Phase 1 — M1: the validated off-target core (must ship — this is the demo)**
- [ ] Reference set: SMILES + ChEMBL off-targets + DILIrank/DICTrank organ join — start with the auto-derivable off-target/cardiac-CNS slice (7/20 clean per Phase-0); load the **7 verified pairs** as fixtures + a documented **novelty search log** (§2).
- [ ] R4 module across the rich panel (hERG, 5-HT2B + aminergic GPCRs, transporters, enzymes): class-membership scoring (pChEMBL≥6, top-5 Tanimoto, z vs background), leave-one-out.
- [ ] **Evidence-tier + coverage flags as first-class output:** *mechanism-linked / weak-coverage / **abstain***.
- [ ] organ/mechanism → assay mapping → ranked panel + evidence trail; **retrospective leave-one-out demo** (§6).

**Phase 2 — M2: breadth for the metabolism/phenotypic gap**
- [ ] Wrap outcome-models for what R4 can't reach: **DILI QSAR** (DILIrank), **mitochondrial** (Tox21 mito-membrane-potential model), **reactive-metabolite structural alerts**; optionally fold in **ADMET-AI / ProTox** endpoints.
- [ ] **Fuse M1 + M2** into one cross-domain assay panel, each item tagged by evidence tier + confidence.
- [ ] Present liver/mito explicitly as *model-predicted* (not mechanism-linked) — honest tiering, not overclaim.

**Phase 3 — Rigor & hardening (post-MVP)**
- [ ] Measure **per-target AUC beyond hERG** (validate the breadth claim); swap class-Tanimoto → trained target-prediction (SEA / per-target QSAR); metabolite-active handling (fenfluramine).
- [ ] Scaffold-split CV, calibration, **assay-recovery top-1/top-3**, baselines (ECFP, ProTox/ADMET-AI), ablation.
- [ ] Scale curation to ~200 drugs (~2 person-weeks); build the matched No-concern comparator.

**Phase 4 — Stretch**
- [ ] Learned-embedding rung (only with demonstrated lift); off-target docking panel; AACT `why_stopped` mining; 3D O3A/ROCS test.

---

## 6. Demo / pitch (investor-facing)

Claim: **"our AI reasons like an experienced med chemist — links a candidate to clinical failures by *shared mechanism* (not just lookalikes), shows why, recommends the next experiment, and we can show it generalizes."** Framed honestly as a **retrospective leave-one-out demo** — hide a known failure's label, query it, show the system recovers the liability. Built on a **verified** pair (Terfenadine × Cisapride):

1. Paste **terfenadine** as the "candidate," its own tox label hidden (leave-one-out).
2. *Screening against N toxicity-causal clinical failures vs. matched No-concern drugs…*
3. **The honest hook:** *"A 2D-similarity tool scores terfenadine against cisapride at 0.19 — it would miss the link. But both block hERG."*
4. **Mechanism result (real, from our prototype):** *"terfenadine flags the hERG class at **z = +6.4** above background — even with terfenadine and cisapride removed from the reference set. Shared predicted hERG liability with drugs withdrawn for QT/TdP."* (An organ-level `N×` enrichment can be shown too — but only once computed against the matched comparator; don't display a placeholder.)
5. **Ranked assay panel** — hERG patch-clamp → iPSC-cardiomyocyte, each tagged with the driving failed drug + tier + citation.
6. **Trust signal:** paste an out-of-domain molecule → the tool **abstains** ("outside applicability domain").
7. Mechanism-edge network; click a node to inspect a failed drug.

**Value prop:** broad safety panels are expensive — *"move the hERG + cardiomyocyte assays to the front of the queue"* saves cost, time, animals, and catches liabilities earlier. Credible because it's grounded in real withdrawals + class-conditional + scaffold-split evidence, not a black box.

---

## 7. Risks & honest framing (raise these before a judge does)

1. **T1 label provenance.** Withdrawn/discontinued ≠ toxicity-caused. → mandatory `failure_reason` + `safety_related` curation before T1 (§4 schema).
2. **Small T1 positive class** (~90–110 in ClinTox). → lean on T2 (DILIrank/DICTrank ~1.3k) for scoring; T1 as the clinical-grounding badge.
3. **"Safe" comparator is a misnomer.** → matched **No-concern** class, not "any approved" (§3e).
4. **Score semantics.** Present enrichment (`N×`) or 0–100 priority, **not** probability-looking decimals unless calibrated.
5. **2D structure is insufficient** (our own test). → mechanism/off-target layer is the engine; state it as a designed choice, not a gap.
6. **Off-target prediction is itself a model** with error — **worst for novel scaffolds** (R4's sensitivity comes from the target's known ligands; measured AUC 0.667 on a diverse isolated-chemotype holdout vs 0.913 pooled). Also: **only hERG has a measured discrimination AUC (0.894)** — the other 36 rich targets are expected-but-unvalidated (Phase 3). → report confidence + applicability-domain flag; measured off-targets for reference drugs, predicted only for the candidate.
7. **Some clinical killers are NOT off-target binding** — reactive/CYP metabolites, idiosyncratic mitochondrial tox, dose/exposure. R4 (M1) physically can't predict these → they are the **M2 outcome-model modules** (DILI QSAR, Tox21-mito, alerts), served at a **lower evidence tier (model-predicted, not mechanism-linked)** and clearly flagged as such. Liver/mito is a *designed module boundary*, not a hidden gap — but M2 is less validated than M1; don't oversell it.
8. **Generalization.** Never quote a random-split number; report scaffold-split CV + AD abstention.
9. **Proprietary gap.** Cleanest discontinued-for-safety data (Pharmaprojects/Citeline, Pharmapendium) is paywalled — note as gold standard if licensed.
10. **Licenses:** DrugBank / SIDER / OFFSIDES non-commercial — flag for any commercialization.

---

## 8. Tools
**Core:** RDKit (ECFP/FCFP, Gobbi pharmacophore, **USRCAT 3D**, conformers via ETKDGv3, descriptors, substructure), **FPSim2** (2D retrieval engine), scikit-learn (scaffold-split CV, calibration, stacking), an LLM for the agent report.
**3D similarity:** RDKit USRCAT / O3A; **Shape-it / Align-it** (open-source ROCS-like); ROCS/Phase (commercial, gold standard).
**Off-target / mechanism (R4):** **Bowes et al. 2012 ~44-target safety panel** (aligned to SAFETYscan47) as the **target vocabulary**; SEA (Similarity Ensemble Approach), ChEMBL target-prediction models, PIDGIN; **per-target QSAR (hERG has strong public models)**; optional docking (hERG/5-HT2B/BSEP panel); reference off-targets = measured ChEMBL bioactivity.
**Baselines / cross-checks:** plain ECFP, ProTox 3.0, ADMET-AI, admetSAR/ADMETlab 3.0, EPA GenRA.
**Only-if-lift-shown:** MolFormer/ChemBERTa/Uni-Mol as an added rung, gated by the scaffold-split holdout.

---

## 9. Key sources
- ClinTox — MoleculeNet arXiv 1703.00564; TDC `single_pred/Tox`
- WITHDRAWN 2.0 — NAR 2024;52:D1503 (v1: NAR 2016;44:D1080)
- Probes & Drugs — probes-drugs.org/compoundsets
- AACT — aact.ctti-clinicaltrials.org; terminated-trials PMC4444136
- DILIrank / DILIst — FDA LTKB; DILIst ScienceDirect S1359644619303824
- DICTrank — FDA (2023); Drug Discov Today DOI 10.1016/j.drudis.2023.103770
- Tox21 / ToxCast — tox21.gov; EPA figshare invitrodb
- SIDER / OFFSIDES / TWOSIDES — nsides.io; Tatonetti Lab
- Closest prior art — Albert, Skinner et al., *"Comprehensive Analysis of Clinically Discontinued Compounds Using an In Vitro Secondary Pharmacology Panel…"*, *ACS Pharmacol Transl Sci* 2025 — **DOI 10.1021/acsptsci.5c00452** (SAFETYscan47; 52 safety-discontinued compounds 2001–2021)
- Secondary-pharmacology safety panel (R4 target vocabulary) — Bowes et al., *Nat Rev Drug Discov* 2012;11:909
- Read-across — EPA GenRA · Predictors — ADMET-AI; ProTox 3.0 (NAR 2024 W513); ADMETlab 3.0 (NAR 2024 W422)

**Verified demo pairs — citations:**
- hERG shared determinants (terfenadine + cisapride) — Saxena et al., PMC2845965
- 5-HT2B valvulopathy shared assay — Setola et al., *Mol Pharmacol* 2003 (PMID 12761331); fenfluramine — Connolly *NEJM* 1997 (PMID 9271479); pergolide — Zanettini *NEJM* 2007 (PMID 17202454); benfluorex/Mediator — Tribouilloy 2011 (PMID 21981882)
- Dual BSEP+mito DILI (troglitazone + nefazodone) — Aleo et al., *Hepatology* 2014;60:1015
- Sitaxentan hepatotox withdrawal — EMA (2010); FIAU trial — IOM review NBK232082; FIAU pol-γ PMID 8622980; perhexiline CPT PMID 8694852
- SMILES/IDs from PubChem & ChEMBL; ECFP4/FCFP Tanimoto computed locally (RDKit, Morgan r=2, 2048-bit)
