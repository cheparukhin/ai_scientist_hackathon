# De-risk probe: metabolite-aware scoring for the fenfluramine 5-HT2B miss

Method: exact reuse of experiments/score.py scoring — mean top-5 Tanimoto (ECFP4,
Morgan r=2, 2048-bit) to a target's ChEMBL active class, as a z vs the same 24-drug
BACKGROUND. Leave-one-out removes the query InChIKey + named mechanistic partners
(fenfluramine/norfenfluramine/benfluorex + pergolide/cabergoline for 5-HT2B;
terfenadine/fexofenadine + cisapride/astemizole for hERG). ChEMBL actives
(pChEMBL>=6, IC50/Ki/Kd): 5-HT2B 1,214 · hERG 1,483 · SERT 1,563.
Background anchors: 5-HT2B mean=0.233 sd=0.058; hERG mean=0.226 sd=0.055.

## Result 1 - fenfluramine parent vs metabolite vs ensemble (z on 5-HT2B)
| species | mean-top5 | z | ~pctile |
|---|---|---|---|
| parent (fenfluramine) | 0.253 | +0.35 | ~64th (reproduces the miss) |
| metabolite (norfenfluramine) | 0.303 | +1.22 | ~89th |
| ensemble = max(parent,metab) | - | +1.22 | ~89th |
| ensemble = mean(parent,metab) | - | +0.78 | ~78th |

Metabolite scoring HELPS but does not fully rescue: +0.35 -> +1.22 (+0.87 / ~3.5x),
still short of a confident hit (terfenadine->hERG = +6.4 for scale). Residual gap =
ECFP small-molecule blind spot (norfenfluramine too small/dissimilar to the ergoline-
heavy 5-HT2B class). self_in_class=False for norfenfluramine -> not a self-match.
Parent ranks hERG (+1.72) ABOVE its true 5-HT2B (+0.35); metabolite equalizes
5-HT2B (+1.22) with hERG (+1.25), stopping the true target from being lowest.

## Result 2 - generalization (one case each direction)
| case | toxic species | parent z | metab z | ens-max | verdict |
|---|---|---|---|---|---|
| fenfluramine -> norfenfluramine (5-HT2B) | metabolite | +0.35 | +1.22 | +1.22 | HELPS |
| benfluorex -> norfenfluramine (5-HT2B) | metabolite | +0.64 | +1.22 | +1.22 | HELPS |
| terfenadine -> fexofenadine (hERG) | PARENT | +6.38 | +5.56 | +6.38 | HURTS / FP trap |

Both 5-HT2B prodrugs (dealkylated to norfenfluramine) improve -> consistent
"metabolite is the toxicophore" signal, not a one-off. Terfenadine counter-example:
parent is the hERG blocker; active metabolite fexofenadine is clinically NON-cardiotoxic
yet still scores +5.56 (near-identical structure) -> metabolite-only would LOWER the
correct signal and ECFP cannot separate detoxified metabolite from toxic parent.
Ensemble-MAX preserves the correct call (max(6.38,5.56)=6.38); MEAN would dilute it.

## Recommended rule
Score the parent AND its known major active metabolite(s); aggregate by MAX (never mean).
MAX lets a toxic metabolite rescue an inactive parent (fenfluramine/benfluorex) while a
detoxified metabolite of a toxic parent cannot dilute the parent's true signal
(terfenadine). Apply metabolite expansion to prodrugs / drugs with a known active
metabolite; under MAX it is neutral-to-helpful and never harmful.

## Honest limits
1. Partial, not full, rescue: +0.35 -> +1.22 is a real ~3.5x lift but still a weak
   positive; fixes wrong-molecule error, not ECFP's small-ligand weakness (SEA/QSAR is
   the deeper fix).
2. Requires external metabolism data (which metabolite, and that it bears the
   toxicophore) - unavailable for a novel candidate without CYP site-of-metabolism
   prediction, which adds its own error. Input-availability limit, not scoring limit.
3. Metabolite scoring can mislead: a detoxified near-analog of a toxic parent
   (fexofenadine) yields a high false-positive score; MAX avoids diluting a true signal
   but does not suppress such a false positive.

## Files
- met_experiment.py  runnable (fetch + score; reuses score.py logic)
- met_actives.json   ChEMBL cache (5-HT2B, hERG, SERT)
- results.json       all z-scores, per-target context, ensembles
