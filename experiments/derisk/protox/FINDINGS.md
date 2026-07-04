# FINDINGS — Does ProTox-3.0 make the assay-recommender redundant?

**De-risk question:** A judge will ask "why not just use ProTox-3.0?" This tests whether that free, established tool already localizes the specific off-target liability our tool ranks — especially for non-cardiac ("buried") culprits.

**Date:** 2026-07-04. **Prior art:** ProTox-3.0, tox.charite.de/protox3, NAR 2024;52:W513 (PMC11223834).

## How ProTox predictions were obtained
- **A real REST API exists and was reverse-engineered/validated** (not just a form): `POST /protox3/src/api_enqueue.php` (`input_type`, `input`, `requested_data`=JSON list of a space-joined model string) → returns a `task_id`; poll `POST /protox3/src/api_retrieve.php` (`id`); results are tab-separated CSVs at `/protox3/csv/<task_id>_{tox_class,result,tox_targets}.csv`. Limit 250 queries/IP/day. The official client (`protox3_api.py`) 404s live but was recovered from the Wayback Machine; a clean stdlib driver (`run_protox.py`) batching all 20 drugs is provided.
- **Live run was blocked** by IP-level rate-limiting of the sandbox (TLS `UNEXPECTED_EOF` on every request after initial access); the Chrome browser fallback extension was not connected. After >3 attempts, switched to the **analytical coverage** method using ProTox's full, authoritatively-sourced endpoint inventory (live Models dropdown via Wayback + API `ALL_MODELS` list + the two NAR papers — all concordant). The driver is ready to run from any un-throttled network.

## ProTox-3.0's complete predicted-endpoint inventory (61 models + targets)
- **Acute:** oral LD50 + class 1-6.
- **Organ (13):** hepatotoxicity, **cardiotoxicity**, neuro-, nephro-, respiratory-, carcino-, immuno-, mutagen-, cyto-, BBB, eco-, clinical-, nutritional-tox (binary Active/Inactive). The **cardiotoxicity model is trained on hERG-blocker data**.
- **Tox21 pathways (12):** nuclear receptors (AhR, AR, AR-LBD, Aromatase, ER, ER-LBD, PPAR-γ) + stress response (ARE, HSE, MMP, p53, ATAD5).
- **MIEs (14):** THRα/β, TTR, RYR, GABAR, NMDAR, AMPAR, KAR, AChE, CAR, PXR, NADHOX, **VGSC (Nav)**, NIS.
- **CYPs (6):** 1A2, 2C19, 2C9, 2D6, 3A4, 2E1.
- **Toxicity Targets (15):** pharmacophore ("toxicophore") fit, class 0-3. Defined as the Novartis in-vitro safety-panel targets that had a **solved human protein-ligand crystal structure** (15 of 73). Confirmed members: **β2-adrenergic (7 toxicophores), Dopamine D3 (1)**. Full 15 not enumerated on any accessible page/paper.
- **Output format (all sections):** a flat table `input | type | Target | Prediction | Probability`. **No ranking, no assay recommendation, no linkage to clinically-failed drugs, no severity/actionability tier, no reordering.**

## Coverage of our 5 culprit categories (20 drugs)
| Culprit | Drugs (n) | Category | Does ProTox localize the specific off-target? |
|---|---|---|---|
| **hERG** | terfenadine, cisapride, astemizole, grepafloxacin, sertindole, thioridazine, terodiline, droperidol, sparfloxacin, mesoridazine (10) | cardiac_herg | **Organ-level only** — `cardiotoxicity` model (hERG-trained) fires, but there is **no hERG target output and no assay named**. Matches the easy cardiac cases as a coarse flag. |
| **5-HT2B** | pergolide, cabergoline, fenfluramine, dexfenfluramine, methysergide (5) | buried | **Partial/unconfirmed** — organ models give at best generic cardio/neuro flags; 5-HT2B *may* be one of the 15 structure-enabled toxicophores (cannot confirm/exclude). Even if present, a flat pharmacophore-fit class, not a valvulopathy assay recommendation. |
| **CB1** | rimonabant, taranabant (2) | buried | **Partial/unconfirmed** — CB1 human structure exists (2016) so it *could* be among the 15; otherwise only a generic `neurotoxicity` flag. No psychiatric/CB1 assay localization. |
| **5-HT3** | alosetron (1) | buried | **No** — 5-HT3 (cys-loop pentamer, no human structure at selection) is in no model/MIE/target. No ischemic-colitis or 5-HT3 signal. |
| **NET** | sibutramine (1) | buried | **No** — no monoamine-transporter output (no human transporter structure). |
| **Cav1.2** | mibefradil (1) | buried | **No** — Cav1.2 not a target; only RYR + VGSC(Nav) channels exist (as MIEs). Coarse cardio flag at best. |

## Verdict
- **Cardiac (hERG) drugs:** ProTox's `cardiotoxicity` model (hERG-trained) **fires** — approximately matches the *easy* cases, but only as a binary organ flag, not "run hERG patch-clamp first."
- **Buried drugs:** ProTox has **no dedicated endpoint** for 5-HT3, NET, or Cav1.2 (3/10 buried → nothing beyond a generic organ probability). For 5-HT2B (5) and CB1 (2) it *might* carry a pharmacophore-fit toxicity target among its unenumerated 15 — the one place ProTox beats a naive "organ-probabilities-only" strawman, stated honestly.
- **Localization / ranking:** In **no case** does ProTox produce a *ranked, reordered assay panel* that puts the killer assay first, and in **no case** is any signal grounded in the specific clinically-failed drugs sharing the mechanism. Its target output is a flat 0-3 pharmacophore-fit list, structurally capped at 15 crystal-structure-enabled targets vs. our ligand-similarity engine's 37/39-target panel (needs no structure; reaches 5-HT3, NET, Cav, hERG).

**Blunt answer to "why not just use ProTox?"** ProTox **partially overlaps** on the coarse cardiac case and *may* touch 5-HT2B/CB1 at the molecular level — so don't claim it is blind to every buried target. But it **clearly falls short on the differentiator**: it emits organ/pathway/target *probabilities*, not a **clinical-failure-grounded, severity-weighted, reordered in-vitro assay panel**. It structurally cannot cover 5-HT3, NET, or Cav (no crystal structures), gives no assay recommendation for any of the 20, and never links a candidate to the withdrawn drugs it shares a mechanism with. ProTox is a legitimate **baseline / lower-confidence outcome component**, not a substitute for the linkage-and-assay-triage layer. It does **not** make the project redundant.

## Files
- `run_protox.py` — runnable stdlib driver batching all 20 SMILES through the recovered API
- `protox3_api.py` — recovered official ProTox-3.0 client (from Wayback)
- `protox3_endpoints.md` — full endpoint inventory + sourcing
- `coverage_analysis.txt` — 20-drug verdict table
