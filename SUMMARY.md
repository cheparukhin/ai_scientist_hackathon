# Toxicity → Assay Recommender — one-page summary

*Start here: [`README.md`](README.md) · full detail: [`toxicity-assay-recommender.md`](toxicity-assay-recommender.md) · evidence: [`experiments/`](experiments/) + [`experiments/derisk/`](experiments/derisk/) · novelty: [`novelty-search-log.md`](novelty-search-log.md)*

## What it is
Paste a drug candidate's structure → the tool **reorders your in-vitro safety-assay plan so the experiment most likely to kill the program runs first.** It does this by linking the candidate to **drugs that actually failed in the clinic for toxicity**, via the off-target mechanism they share — including liabilities no structure search would find. It prioritizes experiments; it does **not** claim a molecule is safe.

## The problem — time and cost to kill
Most candidates die, and the expensive waste is *how long and how much you spend before the kill decision*. A broad safety panel is run in a generic default order; if a candidate is doomed by a specific off-target liability, the assay that proves it is often not the one you run first. We move it to the front → **fewer experiments, less time, less cost, fewer animals, earlier go/no-go.** The product metric is **assays until culprit** (how many experiments until you hit the real liability), and we lower it.

## How it works (3 layers, plain version)
1. **Find the mechanism, not the lookalike.** A candidate rarely *looks* like the drug it will fail alongside (2D similarity of real failure-pairs is near-zero, 0.05–0.19). Instead we use **off-target class matching**: score the candidate against the **whole class of known ligands** for each safety-panel off-target (hERG, 5-HT2B, CB1, 5-HT3, ...). Resembling *some* member of a target's ligand class means the candidate may share liability at that target, even when it resembles no single failed drug.
2. **Weight by consequence, not just likelihood.** A hit isn't automatically a killer. Each off-target carries a **severity / action tag** (no-go / counter-screen / monitor), boosted when the candidate resembles drugs that were *actually withdrawn* for that mechanism.
3. **Recommend the assay + tell the truth about confidence.** Map mechanism → named assay, tagged **mechanism-linked / lower-confidence organ-tox / abstain / known analog**. The validated off-target core covers cardiac/CNS/GI/endocrine/immune liabilities; lower-confidence liver, mitochondrial, and reactive-metabolite modules cover the biology it cannot see directly.

## What the de-risking proved (8 feasibility experiments in `experiments/` + de-risking tracks in `experiments/derisk/`, all on real ChEMBL data)
- **The value is real and quantified — for off-target-mediated failures.** On 20 historical failures, for *buried* off-target liabilities we cut **assays until culprit from ~11.3 → ~3.8** (top-3 recovery **7/10 vs 1/10** for a default panel). **6 genuine non-obvious wins** across three mechanisms: 5-HT2B valvulopathy (pergolide, cabergoline, methysergide), CB1 psychiatric tox (rimonabant, taranabant), 5-HT3 colitis (alosetron).
- **A blind hold-out qualified the scope.** With ChEMBL selecting *both* the drugs and the culprit labels, the magnitude replicates (4.7 vs 11.4) — but drugs withdrawn *specifically for an off-target we can reach* are a **narrow slice**: most withdrawals are off-panel (hepatotoxicity dominates) or on-target pharmacology. **Frame the value as "for off-target-mediated failures," not "all withdrawn drugs."**
- **It discriminates:** ROC-AUC **0.89 (hERG)**, replicated **SERT 0.95, AChE 0.92, MAO-A 0.88** — no "flags everything" failure. Broad: **37/39** panel targets serviceable.
- **The value is *buried* liabilities, explicitly NOT hERG.** Default panels already front-load hERG; there we add nothing and are sometimes worse. Stated plainly.
- **The engine's edge is conditional on novelty.** Against a candidate resembling a known failure, cheap 2D similarity matches us; off-target class matching wins decisively only on **novel chemotypes** (naive 2D -> 0/20 under strict hold-out; off-target class matching -> 12/20). The tool ships a **known analog flag** to be honest about when it adds value.
- **Not redundant with ProTox** — ProTox emits organ/pathway probabilities, can't localize 5-HT3/NET/Cav, and never produces a ranked assay plan grounded in clinical failures.
- **Novelty white space survives** an adversarial search; differentiate from AbbVie OTSA, Albert 2025 (SAFETYscan47), ToxEvaluator.

## Main selling points
- **Grounded in real clinical failures**, not textbook alerts — "you resemble drugs that *died* for this, run this assay first."
- **A number, not a vibe:** assays until culprit **11 → 4** on real historical failures.
- **Honest by construction:** evidence tiers + abstain + known analog flags turn the limits into trust features (and pre-empt the obvious ablation critique).
- **Cheap to build:** retrieval + light reasoning over public data (ChEMBL, DILIrank/DICTrank); no model training needed for v1.

## Demo (retrospective leave-one-out)
Hide a known failure's label, paste it as a candidate, show 2D similarity misses the link, off-target class matching recovers the liability, and the plan reorders. **Primary leave-one-out demo: rimonabant** (CB1, z=+10.8, default rank 15 -> rank 1) — under strict LOO the CB1 signal survives because ChEMBL has many *other* CB1 antagonists resembling it, so the tool re-derives the liability from molecules that never failed. **Pergolide is best shown in normal mode, NOT demo mode**: normally the plan moves 5-HT2B valvulopathy from the default panel's #11 to our **#2**, and the tool *surfaces* 5-HT2B as the headline liability — because that's where pergolide resembles an *actually-withdrawn* drug — even though its raw similarity is higher to its own *therapeutic* D2 target (which leads the plan at #1). Raw similarity ranks the plan; clinical grounding decides which flag to feature — a clean "which signal matters" story (the app states this in-line). Under strict LOO, though, pergolide's 5-HT2B signal collapses (it was carried entirely by its withdrawn ergot cousins, which LOO also hides), so demo-mode pergolide honestly reports it could not recover the mechanism — the app now says exactly this rather than mislabeling the recovered D2 hit as a win. Terfenadine/hERG is kept only as a clean validation anchor (z=+6.4), **not** the headline. Then a trust beat: paste cyclosporine -> **abstains**; paste a close known analog -> **known analog flag**.

## Main risks (and the honest answer)
- **Curation is the real bottleneck** — labeling *why* a drug failed (organ/mechanism/provenance) is ~2 person-weeks for a 200-drug set. → scope to auto-derivable off-target drugs first; MVP uses a small hand-curated set.
- **No value on hERG**, and **blind to metabolite-driven tox** — measured: the only misses were metabolite-active drugs (fenfluramine, sibutramine). → liver/mito/metabolite findings are lower-confidence modules, clearly flagged, not oversold.
- **Engine only beats dumb similarity on novel chemotypes** — frame the claim conditionally; the known analog flag makes it honest.
- **Novel chemotypes degrade** (AUC 0.91 pooled → 0.67 diverse-isolated). → abstain + confidence flags.
- **Licenses:** DrugBank / SIDER / OFFSIDES non-commercial — flag for commercialization.

## Status & verdict
De-risking **done** (8 experiment tracks) and the working demo lives in [`app/`](app/). **Verdict: build it, with honest bounds.** The value is real, quantified (assays until culprit 11 -> 4), and validated on 20 failures against a defensible baseline and an adversarial ablation. It is strongest for buried off-target liabilities in novel candidates; it is not a hERG win and not a metabolite-tox solution. Current build state: MVP + phase-2 polish are implemented; next work is hardening, curation scale-up, and production-quality validation.
