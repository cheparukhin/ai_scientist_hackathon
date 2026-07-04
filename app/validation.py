"""validation.py - the measured, cited validation numbers for the in-app panel (N4).

ALL numbers here are transcribed verbatim from:
  - experiments/derisk/expanded/FINDINGS.md   (assays-to-culprit, top-3, money wins, ranks)
  - toxicity-assay-recommender.md §4          (per-target AUC, scaffold-split, ablation)

Do NOT invent or edit these without re-deriving from those files. A test
(tests/test_app.py) re-parses the FINDINGS.md per-drug table and asserts
BURIED_RANK_PAIRS matches it exactly.
"""

FINDINGS_CITATION = "experiments/derisk/expanded/FINDINGS.md"
SEC4_CITATION = "toxicity-assay-recommender.md §4"

# --- headline validation metrics (FINDINGS.md, buried n=10) ---
MEAN_ASSAYS_TO_CULPRIT_DEFAULT = "11.3"   # Baseline B
MEAN_ASSAYS_TO_CULPRIT_OURS = "3.8"       # Ours
TOP3_OURS = "7/10"
TOP3_DEFAULT = "1/10"
MONEY_WINS = 6                             # genuine non-obvious buried wins

# --- per-target AUC (§4) ---
PER_TARGET_AUC = {"hERG": 0.89, "SERT": 0.95, "AChE": 0.92, "MAO-A": 0.88}

# --- scaffold-split & ablation (§4) ---
SCAFFOLD_POOLED = 0.913
SCAFFOLD_NOVEL_ISOLATED = 0.667
ABLATION_NAIVE = "0/20"                    # naive baseline
ABLATION_R4 = "12/20"                      # R4 backbone

# --- per-drug killer-assay rank: (drug, default-panel rank [B], our rank) for the 10
#     BURIED drugs, transcribed from the FINDINGS.md table. ---
BURIED_RANK_PAIRS = [
    ("pergolide", 11, 1),
    ("cabergoline", 11, 1),
    ("fenfluramine", 11, 13),
    ("dexfenfluramine", 11, 13),
    ("methysergide", 11, 1),
    ("rimonabant", 15, 1),
    ("taranabant", 15, 1),
    ("alosetron", 12, 1),
    ("sibutramine", 13, 5),
    ("mibefradil", 3, 1),
]
