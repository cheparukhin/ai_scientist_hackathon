# Toxicity Risk → Assay Recommender

**A retrieval tool that scores a small-molecule candidate's toxicity risk by linking it to drugs that failed in the clinic (or were withdrawn) for toxicity, and recommends which in-vitro assays to run first.**

> Positioning in one line: *"We connect your candidate to clinically-failed drugs it shares a toxicity **mechanism** with — including ones no 2D-structure search would find — and rank which in-vitro assays would catch the implicated liability first. It prioritizes experiments; it does not declare a molecule safe."*

---

## 0. Project context & status (read first)

- **Event:** a bio × AI hackathon. Deliverable is an **MVP + investor/judge-facing demo** (see §6), not a validated product.
- **How we got here:** chosen over an alternative — a *cross-specificity / immunogenicity risk tool for biologics* (antibody candidates, e.g. cetuximab/α-gal) — because the small-molecule tox tool is **more tractable**: core datasets are downloadable, and the MVP is retrieval + light reasoning. The biologics idea is a fallback.
- **Key pivot (from our own feasibility test — see §4):** we tested whether cheap **2D structural similarity** actually connects mechanistically-related tox failures. **It does not** (7 verified withdrawn-drug pairs, ECFP4 Tanimoto 0.05–0.19). The shared liability is an off-target *binding* phenomenon 2D structure can't see. So the retrieval engine is **mechanism / off-target similarity**, with 2D fingerprints demoted to a fast "obvious-analog" pre-filter.
- **Repo state:** **greenfield.** `/Users/cheparukhin/hackathon/` contains this doc; a Python venv with RDKit exists in scratch for verification. No product code yet.
- **Guiding principle:** ship a working floor first (2D fingerprint retrieval + class-conditional scoring), but the **differentiator is the mechanism/off-target layer** validated against real clinical failures — *not* a black-box structural model.
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
- **Closest paper** — *"Comprehensive Analysis of Clinically Discontinued Compounds Using an In Vitro Secondary Pharmacology Panel"* (PMC12441832, 2024): **52 small molecules discontinued Phase 1–3 for safety** run through a **47-target / 78-assay** panel. It screens **everything** and argues the primary target doesn't predict off-target risk — it does **not** use similarity to *prioritize a subset* of assays. **That is our niche.**

> **Novelty claim — needs a documented search log, don't overstate.** A *preliminary* search found no public tool/paper combining (clinical-failure grounding) × (similarity/mechanism retrieval) × (specific assay recommendation). Before the pitch, record the exact **databases, dates, query strings, and inclusion criteria** so "we found nothing equivalent" is reproducible rather than an assertion.

**Our two defensible novelty claims:**
1. **Reference set grounded in *toxicity-causal* clinical failures** (not assay-positives, not mere market withdrawals).
2. **Output = a specific, ranked assay panel**, made credible by organ+mechanism-resolved labels and a validation protocol.

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

| Pair | Shared mechanism → organ | ECFP4 | FCFP4 | Both withdrawn |
|---|---|---|---|---|
| **Terfenadine × Cisapride** | hERG block (Tyr652/Phe656) → QT/TdP | 0.19 | 0.15 | ✅ 1998 / 2000 |
| Astemizole × Thioridazine | hERG → QT/TdP | 0.17 | 0.25 | ✅ 1999 / 2005 |
| **Pergolide × Fenfluramine** | 5-HT2B agonism → valve fibrosis | 0.11 | 0.14 | ✅ 2007 / 1997 |
| Cabergoline × Benfluorex | 5-HT2B (norfenfluramine) → valvulopathy | 0.11 | 0.15 | ✅ (Mediator scandal) |
| **Troglitazone × Nefazodone** | BSEP + mitochondrial → DILI | 0.13 | — | ✅ both DILIrank-Most |
| Troglitazone × Sitaxentan | BSEP inhibition → cholestatic DILI | 0.11 | — | ✅ 2000 / 2010 |
| Fialuridine × Perhexiline | mitochondrial → hepatotox | 0.05 | — | ✅ (FIAU: 5 trial deaths) |

Feature/pharmacophore 2D fingerprints (FCFP, Gobbi) gave only **1.3–1.5× lift** (absolute sim still ~0.15–0.25); sanity checks passed (true analogs stayed high: pergolide↔cabergoline 0.41). **Conclusion:** the signal linking these is **shared binding target**, which needs a 3D and/or off-target layer — *not* 2D structure. This is the single most important design fact.

*(Open: a 3D shape/pharmacophore test (USRCAT) on the hERG pair is not yet run — it will tell us if rung 3 helps before committing to rung 4.)*

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
| **R3 3D shape + pharmacophore** (USRCAT; O3A/Shape-it, Align-it) | conformer shape + 3D feature layout | RDKit USRCAT (+ Shape-it/Align-it) | **test-gated** — may catch hERG-type 3D pharmacophores 2D misses |
| **R4 off-target / mechanism** | shared (predicted) protein liabilities: hERG, 5-HT2B, BSEP, mito… | target-prediction (SEA, ChEMBL target-prediction, PIDGIN) or a small **docking panel**; matched to reference drugs' ChEMBL off-targets | **the differentiator** — this is what links the verified pairs |

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
- **Calibration:** does 4× enrichment map to a real elevated failure rate? Report a calibration curve; ranked-but-uncalibrated is fine for *prioritization* if not sold as P(harm).
- **Baselines to beat:** plain ECFP (does R3/R4 add signal?) and ≥1 published predictor (ProTox / ADMET-AI) on the same holdout — demonstrate **lift**.
- **Applicability domain (OECD):** flag and **abstain** outside the reference set's chemical/target space.
- **Ablation:** report each rung's contribution — "R4 is what caught this," measured not asserted.
- **Honest bounds:** structure/mechanism proxies miss metabolism, reactive metabolites, dose/exposure; T1 positive class is small (lean on T2 for scoring, T1 as the clinical-grounding badge).

---

## 5. Roadmap (hackathon-scoped, layered by risk)

**Phase 0 — Curate toxicity-*causal* failures (the real work; not just a merge)**
- [ ] Ingest ClinTox, WITHDRAWN 2.0, Probes&Drugs, DILIrank 2.0, DICTrank, a hERG set; resolve names → canonical SMILES (PubChem/ChEMBL); de-dupe.
- [ ] **Curate `failure_reason` + `safety_related` + `organ_system` + `mechanism` + provenance** for each; assign `label_tier`. Filter WITHDRAWN's non-safety reasons out of T1.
- [ ] Pull reference drugs' **known off-targets** from ChEMBL (hERG/5-HT2B/BSEP IC50 etc.) for the mechanism layer.
- [ ] Build the **matched No-concern comparator** (DILIrank/DICTrank No-concern, matched by scaffold/property/route).
- [ ] Load our **7 verified demo pairs** (§4) as fixtures + a documented **novelty search log** (§2).

**Phase 1 — Floor (must ship)**
- [ ] FPSim2 2D retrieval + **class-conditional enrichment** score (vs matched comparator) with applicability-domain flag.
- [ ] organ/mechanism → assay mapping → ranked panel + evidence trail.

**Phase 2 — Mechanism layer + rigor (the differentiator; do before any "wow")**
- [ ] **R3 3D test (USRCAT)** on the verified pairs — keep only if it beats 2D.
- [ ] **R4 off-target/mechanism linkage** — predict candidate off-target profile; match to reference off-targets; this is what should recover the verified pairs.
- [ ] Validation: scaffold-split CV, per-organ enrichment/precision@k, calibration, baselines (ECFP, ProTox/ADMET-AI), ablation.
- [ ] Attribution (shared mechanism + motif/pose) + LLM report + mechanism-edge network.

**Phase 3 — Stretch**
- [ ] Learned embedding rung — only with demonstrated lift.
- [ ] Docking panel for off-targets; AACT `why_stopped` mining on a slice.

---

## 6. Demo / pitch (investor-facing)

Claim: **"our AI reasons like an experienced med chemist — links a candidate to clinical failures by *shared mechanism* (not just lookalikes), shows why, recommends the next experiment, and we can show it generalizes."** Flow, built on a **verified** pair (Terfenadine × Cisapride):

1. Paste a candidate.
2. *Screening against N toxicity-causal clinical failures vs. matched No-concern drugs…*
3. **The honest hook:** *"A 2D-similarity tool scores this against cisapride at 0.19 — it would miss it. But both block hERG."* → shows the 2D limitation up front, then the mechanism link.
4. **Class-conditional result:** *"cardiac-liability enrichment 4.1× vs matched comparator; driven by shared predicted hERG block with terfenadine & cisapride — both withdrawn for QT/TdP."*
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
6. **Off-target prediction is itself a model** with error. → report confidence; use measured ChEMBL off-targets for reference drugs, predicted for the candidate; validate on the verified pairs.
7. **Structure/mechanism proxies miss metabolism / reactive metabolites / dose.** → prioritization, never a safety verdict.
8. **Generalization.** Never quote a random-split number; report scaffold-split CV + AD abstention.
9. **Proprietary gap.** Cleanest discontinued-for-safety data (Pharmaprojects/Citeline, Pharmapendium) is paywalled — note as gold standard if licensed.
10. **Licenses:** DrugBank / SIDER / OFFSIDES non-commercial — flag for any commercialization.

---

## 8. Tools
**Core:** RDKit (ECFP/FCFP, Gobbi pharmacophore, **USRCAT 3D**, conformers via ETKDGv3, descriptors, substructure), **FPSim2** (2D retrieval engine), scikit-learn (scaffold-split CV, calibration, stacking), an LLM for the agent report.
**3D similarity:** RDKit USRCAT / O3A; **Shape-it / Align-it** (open-source ROCS-like); ROCS/Phase (commercial, gold standard).
**Off-target / mechanism (R4):** SEA, ChEMBL target-prediction models, PIDGIN; optional docking (hERG/5-HT2B/BSEP panel); reference off-targets from ChEMBL bioactivity.
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
- Closest prior art — secondary pharmacology of discontinued compounds, PMC12441832 (2024)
- Read-across — EPA GenRA · Predictors — ADMET-AI; ProTox 3.0 (NAR 2024 W513); ADMETlab 3.0 (NAR 2024 W422)

**Verified demo pairs — citations:**
- hERG shared determinants (terfenadine + cisapride) — Saxena et al., PMC2845965
- 5-HT2B valvulopathy shared assay — Setola et al., *Mol Pharmacol* 2003 (PMID 12761331); fenfluramine — Connolly *NEJM* 1997 (PMID 9271479); pergolide — Zanettini *NEJM* 2007 (PMID 17202454); benfluorex/Mediator — Tribouilloy 2011 (PMID 21981882)
- Dual BSEP+mito DILI (troglitazone + nefazodone) — Aleo et al., *Hepatology* 2014;60:1015
- Sitaxentan hepatotox withdrawal — EMA (2010); FIAU trial — IOM review NBK232082; FIAU pol-γ PMID 8622980; perhexiline CPT PMID 8694852
- SMILES/IDs from PubChem & ChEMBL; ECFP4/FCFP Tanimoto computed locally (RDKit, Morgan r=2, 2048-bit)
