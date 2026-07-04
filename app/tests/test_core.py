"""Acceptance tests for the scoring engine (M1) and ranking/plan (M2).

Run:  .venv/bin/python -m pytest app/tests/test_core.py -q
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import core  # noqa: E402
from core import EXAMPLES, score_candidate  # noqa: E402


def _rank_of(result, target_key):
    order = sorted(result["targets"], key=lambda t: result["targets"][t]["z"], reverse=True)
    return order.index(target_key) + 1


# ---------------- M1 ----------------

def test_rimonabant_cb1_top():
    r = score_candidate(EXAMPLES["rimonabant"])
    assert r["status"] == "ok"
    assert r["best_target"] == "CB1_CHEMBL218", r["best_target"]
    assert r["targets"]["CB1_CHEMBL218"]["z"] > 3.0


def test_pergolide_5ht2b_top3():
    r = score_candidate(EXAMPLES["pergolide"])
    assert r["status"] == "ok"
    assert _rank_of(r, "5HT2B_CHEMBL1833") <= 3


def test_cyclosporine_abstains():
    r = score_candidate(EXAMPLES["cyclosporine"])
    assert r["status"] == "abstain"
    assert r["gate"] in ("A", "B")


def test_terfenadine_known_analog_or_herg_top():
    r = score_candidate(EXAMPLES["terfenadine"])
    assert r["status"] == "ok"
    herg_top = r["best_target"] == "hERG_CHEMBL240"
    assert herg_top or r["flags"]["known_analog"]


def test_no_probability_language_in_result():
    # honesty guardrail: engine returns z-scores/ranks, never a P(harm) in [0,1] labelled as probability
    r = score_candidate(EXAMPLES["rimonabant"])
    assert "probability" not in r
    assert "p_toxic" not in r


def test_novel_molecule_in_domain():
    r = score_candidate(EXAMPLES["novel_arylpiperazine"])
    assert r["status"] == "ok"


# ---------------- M2 ----------------

def test_rimonabant_plan_cb1_rank1_delta():
    from core import build_plan
    r = score_candidate(EXAMPLES["rimonabant"])
    plan = build_plan(r)
    cb1 = next(row for row in plan["rows"] if row["target_key"] == "CB1_CHEMBL218")
    assert cb1["our_rank"] == 1, cb1
    assert cb1["default_rank"] >= 12, cb1  # buried by default (~15)
    assert cb1["action"] in ("no-go", "counter-screen")
    assert cb1["evidence"], "CB1 row must carry linked failed-drug evidence"
    # assays-to-culprit headline present
    assert plan["headline"]["our_rank"] == 1
    assert plan["headline"]["default_rank"] >= 12


def test_no_pharm_probability_in_plan():
    from core import build_plan
    r = score_candidate(EXAMPLES["rimonabant"])
    plan = build_plan(r)
    for row in plan["rows"]:
        # priority is a 0-100 integer, never a 0-1 probability
        assert isinstance(row["priority"], int)
        assert 0 <= row["priority"] <= 100
