# Toxicity Risk → Assay Recommender

**A retrieval-based tool that scores a small-molecule candidate's toxicity risk by its similarity to drugs that failed in the clinic (or were withdrawn) for toxicity, and recommends which in-vitro assays to run first.**

> Positioning in one line: *"We retrieve the clinically-failed drugs your candidate most resembles — including non-obvious analogs — show the shared structural motif, and rank which in-vitro assays would catch the implicated liability first. It prioritizes experiments; it does not declare a molecule safe."*

---

## 0. Project context & status (read first)

- **Event:** a bio × AI hackathon. Deliverable is an **MVP + investor/judge-facing demo** (see §6), not a validated product.
- **How we got here:** this idea was chosen over an alternative — a *cross-specificity / immunogenicity risk tool for biologics* (antibody candidates, population-specific pre-sensitization, e.g. cetuximab/α-gal) — because the small-molecule tox tool is **more tractable**: all core datasets are downloadable CSVs and the MVP is retrieval + light reasoning with no model training required. The biologics idea is a possible fallback, not the current plan.
- **Repo state:** **greenfield.** As of writing, `/Users/cheparukhin/hackathon/` contains only this doc — no code yet. Stack/infra not yet chosen.
- **Guiding principle:** ship the **floor** (§4 Layer 1 ECFP4 baseline) first so there is always a working demo; treat learned embeddings, 3D, and AACT mining as differentiators/stretch.
- **Parameters still to fill (TBD — not yet decided in-conversation):**
  - Timeframe / submission deadline.
  - Team size and strengths (eng vs. cheminformatics vs. bio) — affects how far up the stretch layers to reach.
  - Explicit judging criteria (working product vs. scientific novelty vs. commercial potential).
  - Target LLM + vector store + hosting for the demo.
- **Data caveat:** dataset sizes below are **last-published figures** (some flagged "verify"); several sources grow continuously or come from secondary summaries. **Re-confirm exact counts against the primary source before quoting them in a pitch** (especially ClinTox positive-class count, DrugBank withdrawn count, ToxCast endpoint total).

---

## 1. The idea

- **Input:** a small-molecule drug candidate (SMILES; optionally 3D conformer).
- **Reference database:** molecules with known organ-specific toxicity, weighted toward those that **failed clinical trials or were withdrawn *for toxicity*** (not just assay-positive compounds).
- **Output:**
  1. A **per-modality risk vector** — `{liver: 0.71, cardio: 0.20, hERG: 0.10, kidney: 0.05, ...}`.
  2. A **ranked in-vitro assay panel** — e.g. *"prioritize hERG + iPSC-cardiomyocyte (driven by 0.71 similarity to \[drug X, DICTrank-Most-concern]); then hepatocyte spheroid (0.58)."*
  3. An **evidence trail** — the nearest failed drugs, the organ they failed in, and the **shared structural motif / toxicophore** driving the match.

It is a **decision-support / experiment-prioritization** tool, not a binary safety classifier. That framing is both more scientifically honest (structure alone misses metabolism, off-target, and dose-driven toxicity) and easier to validate and adopt.

---

## 2. Why it's novel (white space + prior art to differentiate from)

The two halves of the idea each exist separately; the **combination is unoccupied**:

- **Endpoint predictors** — ProTox 3.0, admetSAR 3.0 / ADMETlab 3.0, ADMET-AI, DeepTox, DEREK. They predict *toxicity probabilities / structural alerts*, **not "which assay to run,"** and are trained on heterogeneous tox data, **not grounded in clinical-trial failures**.
- **Read-across engines** — EPA **GenRA**, OECD QSAR Toolbox, AMBIT. Analog-based tox inference (methodologically the closest cousin), but reference sets are **environmental/regulatory chemicals, not clinically-failed drugs**, and output isn't a prioritized clinical assay panel.
- **Closest single paper** — *"Comprehensive Analysis of Clinically Discontinued Compounds Using an In Vitro Secondary Pharmacology Panel"* (PMC12441832, 2024): took **52 small molecules discontinued in Phase 1–3 for safety** and ran them through a **47-target / 78-assay** panel. Crucially it screens **everything** and argues the primary target does *not* predict off-target risk — it does **not** use structural similarity to *prioritize a subset* of assays. **That is exactly our niche.**

A combined literature query {chemical similarity + clinical-trial failure + toxicity + recommend + in-vitro assay + read-across} returned **0 results** — suggestive that the exact framing is unpublished.

**Our two defensible novelty claims:**
1. **Reference set grounded in *clinical* toxicity failures** (not assay-positives). ClinTox already operationalizes this; WITHDRAWN + Probes&Drugs strengthen it.
2. **Output = a *specific, ranked assay panel*,** made credible by the organ-resolved label layer (DILIrank→liver assays, DICTrank/hERG→cardiac assays, Tox21/ToxCast pathways→named HTS assays).

The genuine white space is the **linkage layer**: candidate → nearest clinically-failed neighbors → the organ/mechanism they implicated → the specific assay(s) that would have caught them.

---

## 3. Datasets

### 3a. Clinical-failure reference set (the novelty; assemble, don't download whole)
No single public dataset cleanly says "these N molecules died in the clinic *because of toxicity*, here are their SMILES." We build it by merging:

| Dataset | Content | Size | Access | Role |
|---|---|---|---|---|
| **ClinTox** (MoleculeNet / TDC) | FDA-approved vs. **failed clinical trials for toxicity** | ~1,491 cpds w/ SMILES; **positives (tox-failures) only ~90–110** | Free — `deepchem.molnet.load_clintox()`, TDC, HF `zpn/clintox` | **Start here.** Its positives were mined from AACT — i.e. our exact idea, already v0. Class-imbalanced. |
| **WITHDRAWN 2.0** (Charité) | Withdrawn/discontinued drugs, **withdrawal reasons classified by organ toxicity** | 578 (v1); 2.0 adds FAERS + mechanism | Free — cheminfo.charite.de/withdrawn; NAR 2024 52:D1503 | **Best organ-resolved failure set.** Higher label quality than ClinTox. |
| **Probes & Drugs** — withdrawn set | Withdrawn drugs pre-joined to ChEMBL/DrugBank structures | ~633 (v02.2025) | Free — probes-drugs.org/compoundsets | De-dupe / canonical-SMILES join across sources. |
| **DrugBank** (withdrawn flag) | `withdrawn` category, full structures | ~250–300 (verify) | Free academic; **non-commercial license** | Structure + flag; weaker on "reason = tox". |
| **AACT** (ClinicalTrials.gov dump) | `why_stopped` free-text, `TERMINATED/SUSPENDED`, interventions | full CT.gov | Free — aact.ctti-clinicaltrials.org | **Stretch goal.** NLP-classify tox terminations → resolve drug names to SMILES (~30–60% resolvable). Demo on a slice; don't rebuild wholesale. |

### 3b. Organ-specific toxicity labels (power the "which assay" recommendation)
FDA-curated *ranked* lists — the highest-quality organ-resolved labels; the backbone of the per-modality score.

| Dataset | Content | Size | Assay it points to |
|---|---|---|---|
| **DILIrank** (FDA LTKB) | Drug-induced **liver injury**, 4 concern classes | ~1,036 (v1); **DILIrank 2.0 = 1,336** | Hepatocyte / spheroid, BSEP inhibition, mito-tox |
| **DILIst** (FDA) | Binary hepatotoxic (broader) | 1,279 | same |
| **DICTrank** (FDA, 2023) | Drug-induced **cardiotoxicity**, 4 concern levels | 1,318 (1,291 approved + 27 withdrawn) | hERG, iPSC-cardiomyocyte, Nav1.5/Cav1.2 |
| **hERG datasets** | hERG channel block | TDC 648; DeepHIT ~14k; UnihERG ~20k | hERG patch-clamp / automated electrophysiology |
| **SIDER** | Side effects mined from labels | 1,556 drugs, 5,868 SE terms | phenotype → organ → assay |
| **OFFSIDES / TWOSIDES** | FAERS-derived off-label / DDI side effects | 438,801 / 868k assocs | pharmacovigilance signal layer |
| **LiverTox** (NIH) | Narrative hepatotoxicity monographs | ~1,200 agents | qualitative liver context (text) |

> ⚠️ SIDER/OFFSIDES/FAERS/DILIrank describe drugs that **reached market** (survivorship bias). Use them as **organ-tox labels**, not as the "clinical-failure" reference set.

### 3c. Assay / mechanism datasets (map similarity → a named assay)
| Dataset | Content | Size | Role |
|---|---|---|---|
| **Tox21** | 12 pathway assays (nuclear-receptor + stress) | ~7,800 cpds × 12 | Each pathway *is* an assay → direct "which assay" signal |
| **ToxCast** (EPA) | High-throughput in-vitro endpoints | ~800+ endpoints × >4,200 chem | Richest public **compound → assay-endpoint** matrix (mechanistic layer). Mostly environmental chemicals — partial drug overlap |
| **ChEMBL** | Measured bioactivities incl. safety targets | ~2.4M cpds | Off-target / safety-pharmacology activities; **MCP available for live pulls** |
| **(from prior notes)** TG-GATEs / Open TG-GATEs / DrugMatrix / LINCS L1000 | Toxicogenomics / histopath / cell response | ~170 / ~600 / >20k | Mechanism enrichment (stretch) |

### 3d. Structure sources (SMILES / fingerprints / name-resolution)
- **PubChem** (public domain, >100M) — best for resolving trial/FDA-list **drug names → structure**.
- **ChEMBL** (CC-BY-SA) — structures + bioactivity; **MCP available**.
- **DrugBank** (free academic, non-commercial) — best annotations.

### 3e. Background / negative set (do not skip)
Approved, chronically-used "safe" drugs (the **ClinTox approved class** works). Required so the score means *"enriched for toxic neighbors relative to safe neighbors,"* not merely *"similar to some molecule."* Without a background, every candidate looks risky because everything has *a* nearest neighbor.

---

## 4. Architecture

Layered so there is always a working floor even if the higher layers slip.

```
             ┌─────────────────────────────────────────────┐
  SMILES  →  │ 1. Similarity channels (interpretable)       │
             │   • ECFP4 Tanimoto (FPSim2)          [floor] │
             │   • pharmacophore fp + physchem distance     │
             │   • 3D shape (USRCAT)             [optional] │
             │   • learned embedding — only if lift proven  │
             └───────────────────┬─────────────────────────┘
                                 ▼
             ┌─────────────────────────────────────────────┐
             │ 2. Retrieval over reference DB               │
             │   • multi-channel k-NN (FPSim2)              │
             │   • neighbors + organ labels + tier + source │
             └───────────────────┬─────────────────────────┘
                                 ▼
             ┌─────────────────────────────────────────────┐
             │ 3. Per-modality scoring (class-conditional)  │
             │   • likelihood ratio: sim(toxic)/sim(safe)   │
             │   • matched background + applicability domain │
             │   → {liver:0.71, cardio:0.20, hERG:0.10,...} │
             └───────────────────┬─────────────────────────┘
                                 ▼
             ┌─────────────────────────────────────────────┐
             │ 4. Assay recommendation (modality→assay map) │
             │   → ranked assay panel w/ driving neighbors  │
             └───────────────────┬─────────────────────────┘
                                 ▼
             ┌─────────────────────────────────────────────┐
             │ 5. Explainability + Agent + Visualization    │
             │   • toxicophore / shared-motif attribution   │
             │   • LLM agent narrates the evidence          │
             │   • clickable toxicity-neighbor network      │
             └─────────────────────────────────────────────┘
```

### Reference DB schema (organ/mechanism-resolved, multi-modal)
```
compound_id | canonical_SMILES | source(s) | provenance/evidence_link | label_tier |
fail_stage | clinical_failure_flag |
liver_risk | cardio_risk | hERG | kidney | neuro | hematologic | ... |
mechanism (optional)
```
- Store labels as a **vector across modalities**, not a single "reason" string — a molecule can carry several liabilities.
- **`label_tier` is non-negotiable for rigor — never silently blend label qualities:**
  - **T1** — withdrawn / discontinued *for toxicity* (WITHDRAWN, ClinTox positives): strongest, but small.
  - **T2** — organ-severity-ranked *marketed* drugs (DILIrank-Most, DICTrank-Most): cleaner labels, but survivorship-biased.
  - **T3** — noisier flags (broad ClinTox-failed, SIDER-derived): supporting evidence only.
- `clinical_failure_flag` (from ClinTox/WITHDRAWN) is a separate badge so we keep the clinical-failure story *and* the score stays credible on the cleaner T2 organ ranks.

### Layer 1 — Similarity channels (interpretable by construction)
Retrieve neighbors with **several complementary, human-inspectable similarity metrics — not one black-box embedding.** Each channel maps to a known toxicity driver, so every retrieval can be explained. Engine: **FPSim2** (fast Tanimoto k-NN over RDKit fingerprints); the reference set is only a few thousand compounds, so all channels run in-memory in real time.

| Channel | Captures | Toxicity rationale | Tool |
|---|---|---|---|
| **ECFP4 / Morgan Tanimoto** *[floor]* | 2D substructure / scaffold | close analogs of known failures | RDKit + FPSim2 |
| **Pharmacophore fingerprint** (ErG / Gobbi 2D-pharmacophore) | feature arrangement (basic amine ↔ aromatic, etc.) | catches **scaffold hops** sharing a pharmacophore (hERG, phospholipidosis) | RDKit |
| **Physicochemical distance** (logP, pKa/charge, PSA, MW) | cationic-amphiphilic / lipophilic character | direct driver of phospholipidosis, hERG, promiscuity | RDKit descriptors |
| **3D shape** (USRCAT) *[optional]* | conformer shape | shape-driven off-target liability | RDKit |

**Why not learned embeddings (the rigorous reason):** an off-the-shelf pretrained encoder (MolFormer/ChemBERTa/Uni-Mol) embeds *general chemistry*, **not toxicity** — there's no guarantee it ranks two hERG blockers closer than fingerprints do, it's unexplainable, and its lift over fingerprints for *this* task is unproven. A true "toxicity space" needs supervised/contrastive fine-tuning on the labels (out of hackathon scope, gated by ~1–3k labeled cpds). **Include a neural channel only if you can demonstrate measurable lift over the interpretable channels on a scaffold-split holdout** (see Rigor). Until that lift is shown, it is unfalsifiable weight — cut it.

### Layer 3 — Per-modality scoring (class-conditional, *not* raw similarity)
**The core rigor point:** similarity to a toxic drug means nothing unless the candidate is *more* similar to toxic drugs than to safe ones — toxic and safe drugs share scaffolds (kinase inhibitors, etc.), so a raw nearest-neighbor score is confounded. Score each organ modality as a **relative enrichment / likelihood ratio** against the *matched* safe background (§3e):

```
LR_m  =  weighted_sim(candidate, toxic-in-m neighbors)
         ─────────────────────────────────────────────
         weighted_sim(candidate, safe-background neighbors)
```

Equivalently: a class-conditional k-NN posterior `P(tox_m | structure)`, or an **enrichment factor** (how enriched the top-k neighbors are for organ-*m* failures vs. the base rate). Combine channels by rank-aggregation or a simple logistic stack — **not** a hidden weighting. Report per organ: the score, the driving neighbors (with provenance + label tier), and an **applicability-domain flag** (§ Rigor). This kills the "everything looks toxic because everything resembles *some* drug" confound.

### Layer 4 — Modality → assay mapping (hand-authored; the payoff)
| Modality flagged | Recommended assays (priority order) |
|---|---|
| Liver (DILIrank) | Human hepatocyte / 3D spheroid → BSEP (bile salt export pump) inhibition → mitochondrial toxicity (Glu/Gal, Seahorse) |
| Cardiac (DICTrank) | hERG patch-clamp → iPSC-cardiomyocyte (contractility / field potential) → Nav1.5 / Cav1.2 |
| Kidney | Proximal tubule cells / kidney organoid → transporter (OAT/OCT) assays |
| Reactive metabolite | Glutathione-trapping / covalent-binding assay; CYP-dependent cytotox |
| Tox21/ToxCast pathway hit | the specific named HTS assay for that pathway (e.g. mitochondrial membrane potential, ARE, ER/AR) |
| Hematologic / neuro / lung | organ-appropriate primary/iPSC model (from neighbor labels) |

### Layer 5 — Explainability, agent, visualization
- **Toxicophore attribution** — return the substructure / motif driving each match (*"shares an aromatic sulfonamide + exposed tertiary amine with cholestatic-injury compounds"*). Methods: fingerprint-bit → atom mapping (trivial for ECFP), or attention / integrated gradients / GNNExplainer for the learned model. **Highest-leverage single feature** — it's the defense against "isn't this just fingerprint similarity?"
- **LLM agent** interprets retrieved evidence into a report (retrieval + reasoning; **no model training needed for v1**). Optional multi-role framing for the pitch: Retriever → Evidence → Hypothesis → Experiment-planner.
- **Visualization** — candidate at center, historical neighbors around it, organ-colored (🔴liver 🔵kidney 🟢heart); click a node → why it failed, shared motif, literature, suggested assay.

### Rigor: validation protocol, calibration & applicability domain
This is what separates the tool from a nice-looking demo. Minimum defensible evaluation, all under **retrospective holdout**:

1. **Scaffold-split cross-validation (never random split).** Random splits leak close analogs between train/test and inflate performance — the tool would look great and generalize to nothing. **Bemis–Murcko scaffold splits** test the honest question: does it flag risk in *novel chemotypes*? Report all metrics under scaffold split.
2. **Per-organ discrimination:** ROC-AUC, and — more relevant for a triage tool — **enrichment factor / precision@k** (of the top-k retrieved neighbors, what fraction truly failed in that organ vs. base rate).
3. **Calibration:** does a 0.7 liver score correspond to a genuinely elevated liver-failure rate? Report a calibration curve. A well-*ranked* but miscalibrated score is acceptable for *prioritization* — as long as you never present it as a probability of harm.
4. **Baselines to beat:** (a) plain ECFP/Tanimoto — does multi-channel actually add signal? and (b) ≥1 published predictor (**ProTox 3.0** / **ADMET-AI**) on the same holdout. A value claim requires demonstrated **lift** over these, not standalone numbers.
5. **Applicability domain (OECD QSAR principle):** flag and **abstain** when a query lies outside the reference set's chemical space (e.g. max Tanimoto to any reference < threshold, or descriptors outside the training range). Rigor = knowing when *not* to answer.
6. **Ablation:** report each similarity channel's individual contribution, so "the pharmacophore channel is what caught this" is measured, not asserted.
7. **Explainability grounded in prior art:** where the shared motif matches a documented **structural alert / toxicophore**, cite it — anchoring the "why" in established SAR rather than model opinion.

**Honest bounds (state these up front):** structure-only → misses metabolism, reactive metabolites, off-target, and dose/exposure-driven toxicity; it is *prioritization*, not a safety verdict; the T1 positive class is small, so lean on T2 for scoring and use T1 as the clinical-grounding badge.

---

## 5. Roadmap (hackathon-scoped, layered by risk)

**Phase 0 — Data assembly (do first; this is the real work)**
- [ ] Ingest ClinTox, WITHDRAWN 2.0, Probes&Drugs withdrawn, DILIrank 2.0, DICTrank, a hERG set.
- [ ] Resolve names → canonical SMILES via PubChem/ChEMBL; de-dupe.
- [ ] Build the modality-labeled reference table with **`label_tier` (T1/T2/T3)** and per-compound provenance.
- [ ] Build a **matched safe background set** — approved, chronically-used drugs, ideally scaffold/property-matched to the toxic set (not just "any approved") so the likelihood ratio isn't confounded by drug-likeness alone.
- [ ] **Curate ≥1 verifiable "non-obvious analog → real clinical tox failure" pair** for the demo (the single most important demo asset; a wrong neighbor when a judge clicks kills the pitch).

**Phase 1 — Floor (must ship)**
- [ ] FPSim2 ECFP4 Tanimoto k-NN retrieval.
- [ ] **Class-conditional per-modality score** (likelihood ratio vs. matched background), with an applicability-domain flag.
- [ ] Modality→assay mapping → ranked assay panel + evidence trail (with neighbor provenance + label tier).

**Phase 2 — Rigor & interpretable channels (this is the defensibility, do before any "wow")**
- [ ] Add pharmacophore + physicochemical similarity channels; rank-aggregate.
- [ ] **Validation:** scaffold-split CV, per-organ ROC-AUC + enrichment/precision@k, calibration curve, applicability domain.
- [ ] **Baselines:** beat plain ECFP and ≥1 published predictor (ProTox / ADMET-AI) on the same holdout; channel ablation.
- [ ] Toxicophore / shared-motif highlighting, cross-referenced to structural alerts.
- [ ] LLM agent report + clickable neighbor network.

**Phase 3 — Stretch (only with demonstrated lift or spare time)**
- [ ] Learned embedding channel — **only if** it beats the interpretable channels on the scaffold-split holdout.
- [ ] USRCAT 3D-shape channel.
- [ ] AACT `why_stopped` tox-termination mining on a slice ("how the reference set scales").

---

## 6. Demo / pitch (investor-facing)

Don't claim "our AI predicts toxicity better than pharma" (needs years of validation). Claim: **"our AI reasons like an experienced med chemist — finds non-obvious precedents, explains them, and recommends the next experiment."** Demo flow:

1. Paste a candidate (e.g. a kinase inhibitor).
2. *Analyzing structure… searching N clinically-failed compounds… found non-obvious cluster.*
3. Reasoning: *"shares a 3D arrangement of aromatic rings + tertiary amine with cholestatic-injury compounds; the scaffold differs, so fingerprint search would miss this."*
4. Ranked assay panel (⭐⭐⭐ hepatocytes → ⭐⭐ BSEP → ⭐ iPSC-cardiomyocyte).
5. Organ-colored neighbor network; click to inspect a failed drug.

**Value prop:** broad safety panels are expensive — *"move renal assays to the front of the queue"* saves cost, time, animal studies, and catches liabilities earlier. Easier to validate than binary safety claims because it **prioritizes tests**.

---

## 7. Risks & honest framing (raise these before a judge does)

1. **Small positive class.** ClinTox has only ~90–110 tox-failures; lean on the FDA organ ranks (DILIrank/DICTrank, ~1.3k each) as the scoring backbone and use clinical-failure flags as a boost.
2. **"Failed in trials" label is noisy** and not always organ-resolved. WITHDRAWN's organ tags mitigate.
3. **Structure alone misses metabolism / reactive metabolites / off-target / dose-driven tox.** → position as *prioritization*, never a safety verdict.
4. **Off-the-shelf embedding ≠ "toxicity space."** Don't over-claim without supervised fine-tuning.
5. **The demo needs a *real* non-obvious analog.** Fabricated clusters collapse under a click-through. Curating a true one is Phase-0 work.
6. **Proprietary gap.** The cleanest "discontinued-for-safety" attrition data (Pharmaprojects/Citeline, Pharmapendium) is paywalled — the 52-compound paper used it. Note as "gold standard if licensed," out of hackathon scope.
7. **DrugBank / SIDER / OFFSIDES licenses are non-commercial** — flag if the output could be commercialized.
8. **Similarity ≠ causation (the confound most tools ignore).** Toxic and safe drugs share scaffolds; a raw nearest-neighbor score is confounded. → addressed by **class-conditional / likelihood-ratio scoring** (§4 Layer 3) against a *matched* safe background.
9. **Generalization to novel scaffolds.** Random-split metrics overstate performance via analog leakage. → **scaffold-split CV + applicability-domain abstention** (§4 Rigor). Never quote a random-split number.

---

## 8. Tools
**Core:** RDKit (Morgan/ECFP + pharmacophore fingerprints + USRCAT shape + physchem descriptors + substructure matching), **FPSim2** (fast Tanimoto k-NN retrieval engine), scikit-learn (scaffold-split CV, calibration, logistic stacking), an LLM for the agent report. **Baselines / cross-checks:** plain ECFP, ProTox 3.0, ADMET-AI, admetSAR/ADMETlab 3.0, EPA GenRA. **Only-if-lift-shown:** HuggingFace MolFormer/ChemBERTa or Uni-Mol as an *additional* channel, gated by the scaffold-split holdout.

---

## 9. Key sources
- ClinTox — MoleculeNet, arXiv 1703.00564; TDC `single_pred/Tox`
- WITHDRAWN 2.0 — NAR 2024;52:D1503 (v1: NAR 2016;44:D1080)
- Probes & Drugs — probes-drugs.org/compoundsets
- AACT — aact.ctti-clinicaltrials.org; terminated-trials analysis PMC4444136
- DILIrank 2.0 / DILIst — FDA LTKB; DILIst ScienceDirect S1359644619303824
- DICTrank — FDA (2023); Drug Discov Today S1359644623002866
- Tox21 / ToxCast — tox21.gov; EPA figshare invitrodb
- SIDER / OFFSIDES / TWOSIDES — nsides.io; Tatonetti Lab
- Closest prior art — secondary-pharmacology of discontinued compounds, PMC12441832 (2024)
- Read-across — EPA GenRA
- Predictors — ADMET-AI (admet.ai.greenstonebio.com); ProTox 3.0 NAR 2024 W513; admetSAR 3.0 / ADMETlab 3.0 NAR 2024 W422
- Encoders — Uni-Mol; MolFormer; ChemBERTa
