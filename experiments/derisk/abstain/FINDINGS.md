# De-risk probe: applicability-domain / "abstain" rule

Goal: make the demo's "paste an out-of-domain molecule → the tool abstains" moment safe
to show live, backed by a concrete, calibrated rule — not a hand-wave.

Data: 7,093 ChEMBL actives across hERG (CHEMBL240), 5-HT2B (CHEMBL1833),
SERT (CHEMBL228), D2 (CHEMBL217), AChE (CHEMBL220). See `ad_actives.json`.

## Key negative finding — the intuitive rule FAILS
The obvious nearest-neighbour rule (abstain if max ECFP4 Tanimoto to any panel active
< threshold) does **not** work with this reference set:
- In-domain drugs span maxNN **0.17–1.00** (metformin 0.18, simvastatin 0.17, prednisone 0.20).
- OOD molecules span maxNN **0.14–0.32** — heavily overlapping the low end of in-domain.
- Any threshold that rejects OOD peptides/sugars would false-abstain ~15 legit drugs.

Root cause: with a narrow 5-target reference set, "low NN similarity" is ambiguous — it
means the same thing for a legitimately low-liability drug and a genuinely out-of-scope
molecule. The best-z score-floor has the identical defect. **Both are demoted to
informational signals.** A low best-z is reported as `weak_coverage` ("no actionable
signal"), which is honestly distinct from "out of applicability domain."

## The rule adopted (abstain if EITHER gate fires)
- **Gate A — scope:** molecule parses AND has ≥6 heavy atoms AND contains no metal.
- **Gate B — OECD descriptor box:** abstain if MW > 800, OR fraction-halogen > 0.30,
  OR fraction-oxygen > 0.33. Thresholds calibrated from reference p99.5, relaxed to the
  observed in-domain/OOD gap.

## Test result — clean on the hand-built set
- **False-abstain: 0/31** — all 6 verified query drugs + 24 background drugs accepted.
- **False-accept: 0/11** — all OOD abstained: cyclosporine, leuprolide (peptide),
  glucose, sucrose, PEG, PFOA, hexachlorobenzene, ethanol, urea, CCl4, cisplatin.

## Honesty caveats
1. Two thresholds are thin: **MW** (levothyroxine 777 in vs leuprolide 816 out, ~5%
   margin) and **fraction-oxygen** (aspirin 0.31 vs PEG 0.36, ~16%). Only the halogen
   and Gate-A margins are wide.
2. This is an "outside-*this-panel*-scope" detector, not a universal drug/non-drug
   oracle — a legit glycopeptide (e.g. vancomycin) would be abstained.
3. The OOD set is hand-built (11 molecules), so 0/0 means "no obvious failure," not a
   guaranteed generalization rate.

## Verdict — safe to demo with a wide-margin molecule
- **Best demo pick: cyclosporine** — a real approved drug but a 1.2 kDa macrocyclic
  peptide (MW 1203 ≫ 800): *"a real drug, but outside our small-molecule off-target
  space, so the tool abstains rather than guesses."*
- Also safe: a sugar, PFOA, or cisplatin.
- **Avoid** the leuprolide fragment (MW 816, barely over the line).
- Do **NOT** describe the abstain as nearest-neighbour-based — that rule was tested and
  does not work here.

## Files
- fetch_ad_actives.py   ChEMBL fetch
- descriptor_probe.py   in-domain vs OOD descriptor distributions
- calibrate_ad.py       threshold calibration + accept/abstain test
- ad_actives.json       ChEMBL cache (5 targets, 7,093 actives)
- ref_pctiles.json      reference descriptor percentiles
- results.json          per-molecule accept/abstain decisions
