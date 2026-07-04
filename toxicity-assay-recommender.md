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
  SMILES  →  │ 1. Representation                            │
             │   • ECFP4 (Morgan r=2) fingerprint  [floor]  │
             │   • pretrained embedding (MolFormer/ChemBERTa)│
             │   • Uni-Mol 3D embedding            [stretch] │
             └───────────────────┬─────────────────────────┘
                                 ▼
             ┌─────────────────────────────────────────────┐
             │ 2. Retrieval over reference DB               │
             │   • Tanimoto / cosine kNN                    │
             │   • returns neighbors + organ labels + source│
             └───────────────────┬─────────────────────────┘
                                 ▼
             ┌─────────────────────────────────────────────┐
             │ 3. Per-modality risk scoring                 │
             │   • similarity-weighted kNN vote per organ   │
             │   • calibrated vs. safe background set       │
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
compound_id | canonical_SMILES | source(s) | fail_stage | clinical_failure_flag |
liver_risk | cardio_risk | hERG | kidney | neuro | hematologic | ... |
mechanism (optional) | representative_evidence_link
```
- Store labels as a **vector across modalities**, not a single "reason" string — a molecule can carry several liabilities.
- `clinical_failure_flag` (from ClinTox/WITHDRAWN) is a separate boost/badge so we keep the clinical-failure story *and* the score stays credible on the cleaner FDA organ ranks.

### Layer 1 — Representation (the "non-obvious similarity" differentiator, honestly scoped)
- **Floor:** ECFP4 + Tanimoto. Guaranteed to run; explainable; ships.
- **Differentiator:** add **one** pretrained encoder — **MolFormer or ChemBERTa** (SMILES, via HuggingFace, no conformers = easiest). Show a case where the learned embedding and fingerprints **disagree** → the "finds what chemists miss" moment.
- **Honesty caveat:** an off-the-shelf pretrained embedding is a *general-chemistry* space, **not a "toxicity space."** True toxicity clustering requires supervised/contrastive fine-tuning on the labels (>1-weekend, gated by ~1–3k labeled cpds). Do **not** over-claim "learned toxicity space" for the MVP.
- **Stretch:** Uni-Mol 3D — only if someone owns conformer generation (3D similarity is conformer-dependent and a time sink).

### Layer 3 — Scoring
Per modality *m*: `score_m = Σ_i sim(candidate, neighbor_i) · label_{i,m}` over top-k neighbors, normalized against the safe-background neighbor distribution (enrichment, not raw similarity).

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

---

## 5. Roadmap (hackathon-scoped, layered by risk)

**Phase 0 — Data assembly (do first; this is the real work)**
- [ ] Ingest ClinTox, WITHDRAWN 2.0, Probes&Drugs withdrawn, DILIrank 2.0, DICTrank, a hERG set.
- [ ] Resolve names → canonical SMILES via PubChem/ChEMBL; de-dupe.
- [ ] Build the modality-labeled reference table + the safe background set (ClinTox approved).
- [ ] **Curate ≥1 verifiable "non-obvious analog → real clinical tox failure" pair** for the demo (see Risks — this is the single most important demo asset; a wrong neighbor when a judge clicks kills the pitch).

**Phase 1 — Floor (must ship)**
- [ ] ECFP4 + Tanimoto kNN retrieval.
- [ ] Per-modality risk vector, background-calibrated.
- [ ] Modality→assay mapping → ranked assay panel + evidence trail.

**Phase 2 — Differentiator**
- [ ] Add MolFormer/ChemBERTa embedding + retrieval; surface a fingerprint-vs-learned disagreement.
- [ ] Toxicophore / shared-motif highlighting.
- [ ] LLM agent report + clickable neighbor network.

**Phase 3 — Stretch (only if time)**
- [ ] AACT `why_stopped` tox-termination mining on a slice ("how the reference set scales").
- [ ] Uni-Mol 3D embedding.
- [ ] Cross-check overlay vs. ProTox / ADMET-AI as an independent second opinion.

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

---

## 8. Tools
RDKit (fingerprints + Tanimoto + substructure), HuggingFace (MolFormer/ChemBERTa), Uni-Mol (3D, stretch), a vector store (FAISS or in-memory for ~5k cpds), an LLM for the agent report. Cross-check baselines: ADMET-AI, ProTox 3.0, admetSAR/ADMETlab 3.0, EPA GenRA.

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
