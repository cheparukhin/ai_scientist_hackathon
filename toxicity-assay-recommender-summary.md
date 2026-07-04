# Toxicity Assay Recommender - Short Plan Summary

## One-line idea

Rank the in-vitro safety assays a team should run first for a small-molecule candidate, using evidence from known clinical failures, off-target liability signals, outcome models, and coverage-aware confidence tags.

## What problem it solves

Drug teams often face a broad safety panel and limited time. The tool does not claim to prove a molecule is toxic or safe. It helps answer a more practical question:

**Which assays are most likely to produce an early decision signal: kill, redesign, counter-screen, or monitor?**

## Core plan

The system is a tiered evidence portfolio:

1. **M1: Off-target mechanism linkage**
   - Uses a secondary-pharmacology safety panel.
   - Scores whether the candidate resembles known active ligands for targets such as hERG, 5-HT2B, SERT, AChE, MAO-A, and other rich safety targets.
   - Links predicted target liabilities to drugs that failed or were withdrawn through the same mechanism.
   - Output tier: **mechanism-linked**.

2. **M2: Outcome and alert models**
   - Covers areas where off-target similarity is weak: DILI, mitochondrial toxicity, reactive metabolites, and phenotypic pathway assays.
   - Uses DILIrank/Tox21/ToxCast-style models, structural alerts, and optional ADMET-AI/ProTox endpoints.
   - Output tier: **model-predicted / alert-based**.

3. **Fusion and ranking**
   - Combines evidence strength, target severity, clinical-failure grounding, exposure-margin context when available, assay actionability, and coverage flags.
   - Produces a ranked assay panel with labels such as **no-go**, **counter-screen**, **monitor**, **weak coverage**, or **abstain**.

## Main selling points

- **It improves experiment planning, not just tox scoring.**
  The headline metric is whether the tool moves the assay that would catch a fatal liability into the top 1-3 tests.

- **It is grounded in real clinical failures.**
  The reference layer prioritizes toxicity-causal failures and withdrawals, not only assay-positive compounds.

- **It is broad where public data is actually broad.**
  The off-target engine covers much of the secondary-pharmacology safety panel across cardiac, CNS, GI, endocrine, immune, and respiratory liabilities.

- **It is honest about weaker areas.**
  Liver, mitochondrial, reactive-metabolite, dose/exposure, and idiosyncratic risks are lower-confidence modules, not overclaimed mechanism links.

- **It explains recommendations.**
  Each assay is tied to a mechanism, failed-drug precedent, evidence tier, confidence level, and actionability label.

## Demo angle

Run a retrospective leave-one-out demo:

1. Hide the known failure label for a historical failed drug.
2. Paste it as a candidate.
3. Show that 2D similarity alone misses the clinically relevant relationship.
4. Show the mechanism/off-target engine recovers the liability.
5. Show the ranked assay plan.
6. Report the key value metric:

**The killer assay moved from rank N in a naive/default plan to top 3 in our plan.**

Example: terfenadine should move hERG patch-clamp and iPSC-cardiomyocyte assays to the front because the candidate lands in hERG-liability space linked to withdrawn QT/TdP drugs.

## Key validation metrics

- **Top-1 / top-3 killer-assay recovery** on blinded historical failures.
- **Mean assay-rank improvement** versus a naive broad panel or generic safety workflow.
- **Time/cost to first decision signal** compared with an unranked assay plan.
- Supporting metrics: per-target AUC, scaffold-split performance, precision@k, calibration, ablation, and abstention rate.

## Main risks

- **Curation bottleneck:** toxicity-causal failure reason, organ, mechanism, evidence link, and label tier need manual or semi-manual curation.
- **M2 is weaker than M1:** DILI, mitochondrial, reactive-metabolite, and idiosyncratic risks are less directly mechanistic and should be tagged lower confidence.
- **Novel chemotypes:** off-target similarity can degrade when a molecule is far from known ligand space, so applicability-domain flags and abstention are essential.
- **Severity scoring needs expert review:** the no-go / counter-screen / monitor table must be credible, not arbitrary.
- **Exposure margin is often missing:** off-target potency only becomes a kill signal relative to efficacious exposure.
- **Licensing and provenance:** DrugBank, SIDER, and OFFSIDES may constrain commercial use; withdrawn status must not be confused with toxicity-causal failure.

## Soundness score

**8.2 / 10**

The plan is scientifically plausible and product-relevant if positioned as assay triage under uncertainty. The biggest thing to prove is not absolute toxicity prediction, but that the system consistently ranks decision-relevant assays earlier than a default plan.
