# How we arrived at this method

*The design story in one page. Every number traces to a runnable experiment in
[`experiments/`](experiments/README.md); the pitch is in [`SUMMARY.md`](SUMMARY.md).*

## The problem
Most drug candidates die. The expensive waste isn't that they die — it's **how long you screen
before you find the reason.** A broad in-vitro safety panel is run in a generic default order; if a
candidate is doomed by one specific off-target, the assay that proves it is often *not* the one you
run first. Our job is a ranking: **move the program-ending test to the front.** The metric is
*assays until culprit* — how many experiments until you hit the real liability — and we lower it.

## The obvious idea, and why it fails
The intuitive approach: *"does my candidate look like a drug that already failed?"* — compare it,
structure to structure, against known failures.

**It doesn't work, and we measured why (E1).** A candidate rarely looks like the specific drug it
will fail *alongside*: across 7 verified withdrawn-drug failure-pairs, 2D similarity was only
**0.05–0.19**. Structurally diverse molecules can share the same off-target. Drug-to-drug matching
misses the mechanism.

## Two fancier fixes we tried — and rejected on evidence
Before inventing anything, we tested whether richer similarity would rescue drug-to-drug matching:

- **Pharmacophore / feature fingerprints** (FCFP4, Gobbi; **E2**) — lift of only **1.3–1.5×**,
  absolute similarity still ~0.15–0.25. Doesn't rescue.
- **3D shape matching** (USRCAT, 30 conformers per molecule; **E3**) — improved retrieval rank for
  **4/7** pairs but **worsened 3/7**, one clean rescue. A weak, unreliable supplement at ~30× the
  compute. **Tested and not adopted** — the engine has never used 3D.

The lesson: the problem wasn't the *similarity metric*. It was the *question*.

## What worked: off-target class matching (E4)
Keep the same cheap 2D fingerprint (ECFP4 Tanimoto) — but change what you compare against. Instead
of *candidate → single failed drug*, score **candidate → the whole class of molecules known to hit
each off-target** (its ChEMBL ligands), taking the mean of the top-5 similarities.

A candidate can resemble **no single failed drug** yet still match a target's *ligand class* — and
that class match is what flags the shared mechanism. Validated in E4: comparing to the class,
terfenadine scores **0.57 (z = +6.4)** at hERG versus **0.19** pairwise; across 6 drugs the correct
mechanism was the top-scoring of the 3 tested targets for **4/6** and above background for **all 6**.

**Worked example — rimonabant at CB1** (reproducible in-app under *Show the calculation*):

| | Value |
|---|---|
| Known CB1 binders compared (ChEMBL) | 1,266 |
| Mean top-5 Tanimoto to that class (raw) | 0.940 |
| Background (25 ordinary drugs): mean / SD | 0.229 / 0.062 |
| **z = (0.940 − 0.229) / 0.062** | **+11.4** |
| Default panel rank → ours | **#15 → #1** |

The proof it's mechanism, not memorization: hide rimonabant *and its withdrawn cousins* (Demo mode),
and the nearest single failed drug drops to **0.17** — naive matching sees nothing — yet class
matching still recovers CB1 at **z ≈ +10.8** from *other* CB1 ligands that never failed.

## From a score to a decision
Two more layers turn per-target enrichment into an assay plan:

- **Consequence, not just likelihood.** Two signals do different jobs. A per-target **severity
  tier** gently reweights the plan order (a high-severity target isn't buried under a marginally
  higher-scoring low-severity one). Separately, the **headline** — the one liability we feature — is
  the highest-priority target where the candidate resembles a drug that *actually failed in the
  clinic* (grounding), which need not be the plan's #1. Example: pergolide's raw similarity is
  highest to its own *therapeutic* **D2** target, so D2 leads the plan; but the tool **headlines
  5-HT2B valvulopathy** (our #2, vs the default panel's #11) because that's where pergolide
  resembles an actually-withdrawn drug. Raw similarity ranks the plan; grounding chooses the flag.
- **Map to a named assay** and tell the truth about confidence: mechanism-linked vs lower-confidence
  organ-tox modules vs **abstain** (out of applicability domain) vs **known-analog** (a plain lookup
  would already catch it).

## What the validation proved — and its honest bounds
- **The win, quantified (n=20 historical failures):** for buried off-target liabilities, *assays
  until culprit* fell from **11.3 → 3.8**; the killer test landed top-3 for **7/10** candidates vs
  **1/10** for a default panel; **6 non-obvious wins** a default order misses. A blind hold-out
  replicated the magnitude (**4.7 vs 11.4**).
- **It discriminates (not a flag-everything model):** ROC-AUC **0.89 hERG**, replicated **SERT 0.95,
  AChE 0.92, MAO-A 0.88** (E5, E8). Broad reach: **37/39** panel targets serviceable (E7a).
- **Where it's weakest, stated plainly:** the edge is *buried* off-targets, **not hERG** (panels
  already run hERG first); it's **blind to liver and metabolite-driven toxicity** (liver is a data
  desert — E7a); and it degrades on truly novel isolated chemotypes (scaffold-split AUC
  **0.913 → 0.667**, E6). Those cases surface as lower-confidence modules, abstain, or a
  known-analog flag — never a false "clear."

## The one-sentence version
> We kept cheap 2D similarity but pointed it at the right target — a candidate's resemblance to the
> *whole ligand class* of each safety off-target, weighted by clinical consequence — to move the
> program-ending assay from ~#11 to ~#4, honestly, for off-target-mediated failures.
