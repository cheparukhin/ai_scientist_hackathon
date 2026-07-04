# Assay-recovery — THE crown-jewel metric, computed for real

**Question.** "Our tool moves the assay that would have caught the liability from rank N in a
naive/default safety-assay plan into the top-3." Computed honestly on the 6 verified failed
drugs against a 14-target safety panel.

## Method
- **Panel:** 14 ChEMBL safety targets, each → a named assay. Actives at pChEMBL>=6,
  standard_type in {IC50,Ki,Kd}. Counts: hERG 1483, Nav1.5 258, Cav1.2 70, 5-HT2B 1214,
  5-HT2A 1505, D2 1523, M2 1397, alpha1A 1354, mu-opioid 1394, SERT 1563, AChE 1310,
  MAO-A 831, COX-1 393, H1 1285. (ar_actives.json)
- **Scoring:** identical to score.py — per target, mean-top5 ECFP4 (Morgan r=2, 2048-bit)
  Tanimoto to that target's actives, z-scored vs the same 24-drug background. Strict leave-one-out:
  query InChIKey + its mechanistic-class partners removed from every class before scoring.
- **Three rankings of the 14 assays, killer-assay rank recorded under each:**
  - Ours — assays by candidate z-score, descending.
  - Baseline A (alphabetical) — the strawman the doc uses.
  - Baseline B (defensible, candidate-INDEPENDENT):
    hERG -> Nav1.5 -> Cav1.2 -> 5-HT2A -> 5-HT2B -> D2 -> alpha1A -> M2 -> H1 -> mu-opioid ->
    SERT -> MAO-A -> AChE -> COX-1. Cardiac ion channels front-loaded (hERG always #1), then
    core aminergic GPCRs, then transporters/enzymes. 5-HT2B placed at rank 5 among the
    aminergic GPCRs — a CONSERVATIVE (baseline-favouring) placement, so any win over Baseline B
    is understated, not inflated.

## Result — killer-assay rank (of 14)

| Drug | Killer assay | Ours | Baseline A (alpha) | Baseline B (default) | killer z |
|---|---|---:|---:|---:|---:|
| terfenadine | hERG patch-clamp | 1 | 9 | 1 | +6.38 |
| cisapride | hERG patch-clamp | 1 | 9 | 1 | +5.22 |
| astemizole | hERG patch-clamp | 1 | 9 | 1 | +6.22 |
| thioridazine | hERG patch-clamp | 3 | 9 | 1 | +2.59 |
| pergolide | 5-HT2B counter-screen | 1 | 2 | 5 | +2.50 |
| fenfluramine | 5-HT2B counter-screen | 11 | 2 | 5 | +0.35 |

### Recovery rates (n=6)
|  | Top-1 | Top-3 |
|---|---:|---:|
| Ours | 4/6 | 5/6 |
| Baseline A (alphabetical) | 0/6 | 2/6 |
| Baseline B (defensible default) | 4/6 | 4/6 |

## The honest read — does the headline survive a defensible baseline?

Against the alphabetical strawman (Baseline A): yes, easily (hERG 9 -> 1; top-3 5/6 vs 2/6).
But alphabetical has no clinical meaning and, by accident of the "5" in its name, even
front-loads 5-HT2B to rank 2 — so the strawman is not even consistently bad.

Against the defensible baseline (Baseline B): mostly NO.
- The four hERG drugs get nothing from us. A real secondary-pharmacology panel front-loads hERG
  at rank 1 already. For terfenadine/cisapride/astemizole we tie (1=1); for thioridazine we are
  actually WORSE — our tool ranks hERG at 3 (its D2 on-target z=+4.83 and 5-HT2B z=+2.60 outrank
  hERG z=+2.59) while the default panel had hERG at 1. So "we move hERG into the top-3" is an
  EMPTY claim for hERG drugs: it was already there.
- Top-1 recovery is a tie (4/6 = 4/6). Top-3 differs by exactly one drug: pergolide.

The one genuine, non-obvious win — pergolide. 5-HT2B (valvulopathy) is a specialized
counter-screen sitting at rank 5 in the default panel — behind three cardiac ion channels and
5-HT2A. Our tool surfaces it to rank 1, placing the fibrotic-valve killer above pergolide's own
dopaminergic on-target (5-HT2B z=+2.50 > D2 z=+2.14 > 5-HT2A z=+1.63). Default rank 5 -> ours
rank 1: a real move into the top-3, exactly the reordering the pitch needs.

The honest miss — fenfluramine. Killer 5-HT2B lands at rank 11, WORSE than both baselines
(z=+0.35, background-level). Expected and disclosed: the real 5-HT2B binder is the metabolite
norfenfluramine, so the parent carries no signal. Do not demo this one.

## Verdict
Stated generically ("we move the killer assay into the top-3"), the crown-jewel claim DOES NOT
survive a defensible standard-practice baseline for the cardiac/hERG cases — hERG is already
rank 1 by default and for thioridazine we are marginally worse. It survives in exactly ONE of
six cases, and it is the compelling one: pergolide -> 5-HT2B, default rank 5 -> ours rank 1.
Honest framing of the headline is NOT "hERG first" (industry already does that) but:
"for a valvulopathy liability that a standard panel buries below the cardiac channels and the
drug's own on-target, we surface the 5-HT2B counter-screen to the front."

Strongest honest demo: pergolide — the only genuinely non-obvious reordering that beats a
defensible baseline. Keep terfenadine/hERG as the validation anchor (z=+6.4 clean recovery) but
do NOT hang the "we reorder the panel" headline on it: it is a no-op against real practice.

## Files
- ar_fetch.py — ChEMBL fetch (-> ar_actives.json)
- ar_score.py — leave-one-out scoring + three-way ranking (-> results.json)
- results.json — full z-vectors, per-drug ranks, recovery rates, both baseline orders
