# Novelty Search Log — Toxicity Risk → Assay Recommender

**Purpose.** Reproducible evidence for (or against) the core novelty claim in [`toxicity-assay-recommender.md`](toxicity-assay-recommender.md) §2. A judge will ask "did you actually search?" — this is the answer.

**The claim being tested (the "triple combination"):**
> No public tool/paper combines **(A) a reference set grounded in *toxicity-causal* clinical failures** × **(B) similarity / mechanism-based (off-target) retrieval** × **(C) a specific, RANKED in-vitro assay recommendation** (which assay to run first, grounded in the linked clinical failures).

**Search dates:** all searches run **2026-07-04**.
**Databases/tools:** Amass BiomedCore (≈40M PubMed/PMC records), web search (WebSearch, US index) for commercial software/platforms not in the peer-reviewed literature. Amass BiomedCore does **not** support an exact per-query hit count (returns a top-ranked page); "#hits reviewed" below = number of records inspected from the ranked result set. This is the reproducibility limitation to disclose.

**Bottom line:** The triple-combination white space **SURVIVES** the search. No published paper or software tool was found that does all three of A×B×C. However, **every pairwise combination is occupied**, and several tools are "one step away." The honest novelty is a **novel integration**, not a novel capability — pitch it that way.

---

## 1. Search table

| # | Tool / DB | Query string (verbatim) | #reviewed | Verdict |
|---|---|---|---|---|
| 1 | Amass BiomedCore | `clinical trial failure toxicity drug withdrawal machine learning similarity read-across assay recommendation` | 10 | No triple. Returns tox-prediction / ML-withdrawal / review papers (Liu liver-ML; Zatorski withdrawal-ML; AI-tox reviews). None recommends a ranked assay panel. |
| 2 | Amass BiomedCore | `off-target secondary pharmacology safety panel assay prioritization triage compound similarity` | 10 | **Richest vein.** Surfaces the closest prior art on the B axis: Muthas 2013 (pharmacological-similarity → safety), Deaton 2019 (rationalize which targets to screen), ToxEvaluator (Pelletier 2016), OPLE (Biehn 2026), plus Albert 2025. None grounds retrieval in a curated clinical-failure set AND outputs a ranked assay plan. |
| 3 | Amass BiomedCore | `read-across in vitro assay selection recommendation toxicity mechanism` | 10 | Read-across literature is environmental-chemical / REACH / nanomaterial-centric (Zhu 2016, Ball GRAP, EU-ToxRisk, Mizuno 2026 pesticide thyroid). Confirms doc §2: read-across engines exist but reference sets are environmental and output is a hazard call, not a clinical assay panel. |
| 4 | Amass BiomedCore | `recommend prioritize in vitro safety assay panel candidate molecule cardiotoxicity hepatotoxicity ranked` | 10 | Endpoint-prediction platforms (Kadioglu 2021 cardiotox ML + ProTox-II filtering; DILIsym QST) and phenotypic-screening reviews. Predict tox / filter compounds; none *re-orders an assay plan* from clinical-failure linkage. |
| 5 | Amass BiomedCore | `withdrawn drug similarity network link candidate shared mechanism organ toxicity retrieval` | 10 | A-axis hits: Xue 2013 (pharmacological network, withdrawn-vs-approved by shared target → safety indication), Zatorski 2025 (withdrawal ML), Li 2017 (DIM systems biology). Output = safety flag / withdrawal probability, **not** a ranked assay. |
| 6 | Amass BiomedCore | `select which safety assay to run drug candidate decision support experiment prioritization off-target` | 10 | Decision-support systems exist but for other domains (REFINE DSS = nanomedicine ITS → "prioritized lists of assays"; TargetDB tractability). Confirms "prioritized assay list" is a known DSS pattern — but not for small-molecule tox grounded in clinical failures via similarity. |
| 7 | Amass BiomedCore | `connect drug candidate to failed withdrawn drugs shared off-target recommend which assay would have caught liability` | 10 | Xue 2013, Zatorski 2025, Zhao 2012 (withdrawn-CV-drug targets human vs mouse), OPLE, Path4Drug. A×B present in several; C (ranked assay) present in none. |
| 8 | Amass BiomedCore | `get_record DOI 10.1021/acsptsci.5c00452` (Albert/Skinner SAFETYscan47) | 1 | **Verified — see §3.** Confirms it screens ALL compounds through the 47-target panel; does NOT use similarity to prioritize a subset. |
| 9 | WebSearch | `software tool recommend rank in vitro safety assay panel small molecule off-target clinical failure similarity` | 8 | Commercial/tooling landscape: Eurofins SAFETYscan47, Reaction Biology panels, Off-targetP ML, AbbVie OTSA, optimized 50-assay panels. All predict/perform off-target screening; none rank assays from a clinical-failure-grounded similarity link. |
| 10 | WebSearch | `Off-Target Safety Assessment OTSA AbbVie hierarchical similarity clinical alerts predict liabilities small molecule` | 8 | **Closest commercial framework — see §2.** OTSA: hierarchical similarity (2D/SEA/QSAR/3D/docking/ML), >7,000 targets, ~2M preclinical+clinical alerts. Predicts off-target liabilities and links to alerts; not a ranked in-vitro assay *plan*. Also surfaces OFF-X (Clarivate) safety-alert DB. |
| 11 | WebSearch | `retrospective analysis would in vitro assay have detected clinical toxicity failure prioritize secondary pharmacology panel similarity to discontinued drugs` | 8 | Re-confirms Albert 2025 (52 discontinued cpds retrospectively through SAFETYscan47) and the Brennan/Urban preclinical secondary-pharmacology resource (target–ADR associations). Retrospective "would it have caught it" is done by *running the whole panel*, not by similarity-prioritizing a subset. |

---

## 2. Closest prior art — annotated, with precise differentiation

Grouped by which of the three axes each occupies. **The pattern to notice: A×B and B×C both exist; A×B×C→ranked-assay does not.**

### Tier 1 — the must-differentiate items (closest overall)

**1. AbbVie OTSA — "Off-Target Safety Assessment" (in silico off-target profiling framework).**
*What it does:* hierarchical similarity (2D Tanimoto → SEA → local QSAR → 3D pocket similarity → docking → ML) over >7,000 targets, cross-referenced to ~2M preclinical + clinical safety alerts, to predict a candidate's off-target liabilities before/beyond the in-vitro panel.
*How it differs:* This is essentially **our R4 engine, industrialized, plus an alert layer.** But (a) its reference is a broad alert database, **not a curated *toxicity-causal clinical-failure* set** (T1); (b) its output is a **predicted off-target/liability list**, not a **re-ordered in-vitro assay plan** with an actionability tag; (c) it does not frame the deliverable as "the assay that would have caught this moves to top-3." Closest thing to us on A(partial)×B; **misses C**. *Cite it; position our differentiator as the assay-triage + clinical-failure grounding on top of an R4-style engine, which OTSA itself demonstrates is feasible.*

**2. Albert, Skinner et al. 2025 — SAFETYscan47 on 52 discontinued compounds (ACS Pharmacol Transl Sci, DOI 10.1021/acsptsci.5c00452). VERIFIED (§3).**
*What it does:* runs 52 safety-discontinued clinical compounds (2001–2021) through the **full 47-target** Eurofins panel; shows the panel detects off-target activities that rationalize the ADRs, and that the primary target doesn't predict off-target risk.
*How it differs:* occupies **A×C-ish but by brute force** — it *has* the clinical-failure set and *does* map to assays, but it **screens everything** and uses **no similarity/mechanism model to prioritize a subset** of assays for a *new* candidate. It is a retrospective validation of a fixed panel, not a per-candidate recommender. **This gap — similarity-driven assay prioritization — is exactly the project's claimed niche, and the search confirms the gap is real.** Must-cite; it is the doc's own "closest paper" and it holds up.

**3. ToxEvaluator (Pelletier 2016, Pfizer × CTD; Database journal).**
*What it does:* enter a compound + target + adverse-event term; returns structurally similar public + proprietary compounds, known + predicted secondary (off-target) pharmacology, and CTD-curated links between targets and the adverse finding.
*How it differs:* conceptually the nearest **A×B integration** — similarity + off-target + tox-linkage in one tool. But it is a **post-hoc interpretation aid for an already-observed tox finding**, not a prospective recommender; it does **not** produce a **ranked in-vitro assay plan**, and its reference is internal tox studies, not a curated clinical-failure set. Must-differentiate: we are prospective and we output the ranked experiment plan.

### Tier 2 — occupy one or two axes, differentiate briefly

**4. Deaton et al. 2019, "Rationalizing Secondary Pharmacology Screening…" (Tox Sci, DOI 10.1093/toxsci/kfy265).** Uses human genetics + pharmacology to **select which off-target proteins to include in the screen** (i.e. panel prioritization, the C-axis idea). But prioritization is driven by **genetic/phenotype evidence, not similarity to clinically-failed drugs**, and there is no per-candidate assay re-ranking. Cite as prior art for "prioritize the panel."

**5. Muthas & Boyer 2013, "Exploiting Pharmacological Similarity to Identify Safety Concerns" (Mol Inf).** Treats the **secondary-pharmacology profile itself as a similarity descriptor** and reads across the in-vivo findings of "pharmacological neighbors" — explicitly beyond chemical similarity. Strong B-axis precedent and philosophically aligned. But no curated clinical-failure grounding and no ranked assay output.

**6. Zatorski & Schlessinger 2025, "Capturing Unanticipated Drug Toxicities…" (Research Square preprint).** Ensemble ML on protein targets + structure + fingerprints predicts **whether a drug will be withdrawn for toxicity** (92% acc); recovers BSEP as a feature. Occupies **A×B** but output = a **withdrawal probability**, not an assay plan. Closest recent A×B ML paper; cite it.

**7. Xue et al. 2013, "…Comparative Analysis of Withdrawn Drugs Using Pharmacological Network."** Links candidates to **withdrawn drugs by shared target** via a network + 2D fingerprint similarity to give a **safety indication**. A×B again; output is a flag, not an assay. (Note: relies on the very 2D similarity the project's own feasibility test showed fails on mechanistic pairs — a point in the project's favor.)

**8. OPLE (Biehn 2026) / Off-targetP ML / OTSA-class tools.** Predict off-target activity against fixed SafetyScreen panels (18/44/47) from similarity + ML. Pure **B**; no clinical-failure grounding, no assay ranking. These are the *engine* the project would build, not the product — consistent with the doc's framing that "the off-target predictors are the means, not the product."

**9. Read-across engines — EPA GenRA, OECD QSAR Toolbox, AMBIT (confirmed via search #3).** Analog-based tox inference, but reference sets are **environmental chemicals** and output is a hazard/endpoint call, not a clinical in-vitro assay panel. Matches doc §2 exactly.

**10. Endpoint predictors — ProTox 3.0, admetSAR/ADMETlab 3.0, ADMET-AI, DeepTox, DEREK (confirmed as landscape, e.g. Kadioglu 2021 using ProTox-II as a filter).** Output tox probabilities / structural alerts. Not grounded in clinical-trial failures; **do not say "which assay to run."** Matches doc §2.

### Adjacent infrastructure (not competitors — the project's data layer)
- **Bowes et al. 2012** (44-target panel) and the **Brennan/Urban 2023 preclinical secondary-pharmacology resource** (target–ADR associations) — feed the severity/actionability table; not recommenders.
- **OFF-X (Clarivate)**, **Comparative Toxicogenomics Database** — curated safety-alert / target–tox knowledge bases; candidate data sources, not tools that do the triple.

---

## 3. Verification of the Albert 2025 SAFETYscan47 paper

- **DOI 10.1021/acsptsci.5c00452 resolves.** PMID **40969870**; *ACS Pharmacology & Translational Science*, published **2025-09-12**; 15 authors (Albert, Skinner, … Jahic; Eurofins Discovery). Amass record `AMBC_LjubRI0EYQEyLLjKcG0ykl9Uq4c`. Full text on PMC (PMC12441832).
- **What it does (from the verified abstract):** demonstrates the utility of in-vitro secondary-pharmacology profiling by running **clinically failed/discontinued drug candidates** through **SAFETYscan47** (47 targets — GPCRs, ion channels, nuclear receptors, transporters, enzymes) and comparing them against a broader compound-modality dataset, to show the panel surfaces off-target activities that explain the ADRs.
- **Confirmed gap:** it **screens all compounds through the whole panel** and argues the *primary* target does not predict off-target risk. It does **NOT** use chemical/mechanism **similarity to prioritize a subset** of assays for a new candidate, and it does **NOT** produce a per-candidate ranked assay recommendation. The doc's characterization in §2 is **accurate**. This is the single best "closest paper" to cite, and the niche it leaves open is genuine.

---

## 4. What white space survives (and how strong it is)

**Survives:** No public tool/paper does **A (curated toxicity-causal clinical-failure reference) × B (mechanism/off-target similarity retrieval) × C (specific ranked in-vitro assay recommendation, "move the killer assay to top-3")** as a single integrated per-candidate workflow.

**But be honest about the strength (this is the risk):**
- **Every pairwise edge is occupied.** A×B: Xue 2013, Zatorski 2025, ToxEvaluator, OTSA. B×C: Deaton 2019, OTSA, Off-targetP ML, OPLE, REFINE-DSS pattern. A×C (brute force): Albert 2025. The novelty is the **specific triple integration**, not any component — several of which are mature.
- **The assay-recommendation framing is the freshest differentiator.** "Which in-vitro assay to run first, ranked, grounded in the clinical failures it resembles by mechanism" was **not found in any source**. Lead with *that*, plus the assay-recovery metric (rank N → top-3), which no competitor reports.
- **Confidence in the novelty claim: MODERATE-to-HIGH** that no published tool does the exact triple; **MODERATE** that it is defensibly novel to an expert judge, because the integration is only "one step" beyond OTSA + Albert 2025. Frame as **"a novel integration of known components, with a new deliverable (ranked assay plan) and a new grounding (toxicity-causal clinical failures)"** — not "a capability no one has."

## 5. Must-cite / must-differentiate in the pitch

1. **Albert/Skinner 2025 (SAFETYscan47)** — the doc's own closest paper; verified. Differentiate: we *prioritize a subset* via similarity; they screen everything.
2. **AbbVie OTSA** (in-silico off-target profiling, ~2M alerts) — closest *engine*; nearest thing to R4 + clinical alerts. Differentiate: we add clinical-failure grounding + a **ranked assay plan** on top.
3. **ToxEvaluator (Pfizer, 2016)** — closest A×B *integration*. Differentiate: prospective per-candidate recommender vs. their post-hoc finding-interpretation tool.
4. **Deaton 2019** — prior art for "prioritize which secondary-pharmacology targets to screen." Differentiate: our prioritization is similarity-to-clinical-failures, not genetics.
5. **Zatorski 2025 / Xue 2013** — prior art for "link candidate to withdrawn drugs and predict tox." Differentiate: they output a flag/probability; we output the ranked experiment plan.
6. **Read-across (GenRA/QSAR Toolbox) and endpoint predictors (ProTox/ADMET-AI)** — already correctly handled in §2; keep the environmental-chemical / no-assay-output distinctions.

**Reproducibility caveat to state in the pitch:** Amass BiomedCore returns a ranked page, not an exact per-query count, so "#reviewed" is records inspected, not total corpus hits. For a hardened claim, re-run the same query strings against PubMed/Google Scholar/Dimensions with date stamps and record raw hit counts; the qualitative negative result (no A×B×C tool) is expected to hold.
