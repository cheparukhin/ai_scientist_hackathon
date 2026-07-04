"""M2 acceptance (N5): reactive-metabolite alerts + hepatotox/mito read-across.

Run:  .venv/bin/python -m pytest app/tests/test_outcomes.py -q
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from outcome_modules import outcome_panel, outcome_scores, reactive_alerts  # noqa: E402

TROGLITAZONE = "Cc1c(C)c2c(c(C)c1O)CCC(C)(COc1ccc(CC3SC(=O)NC3=O)cc1)O2"
NITROANILINE = "O=[N+]([O-])c1ccc(N)cc1"
METFORMIN = "CN(C)C(=N)N=C(N)N"


def test_troglitazone_hepatotox_and_mito_flagged():
    s = outcome_scores(TROGLITAZONE)
    assert s["hepatotox"]["flagged"] is True, s["hepatotox"]
    assert s["mito"]["flagged"] is True, s["mito"]


def test_nitroaromatic_reactive_alert():
    alerts = reactive_alerts(NITROANILINE)
    assert alerts, "expected reactive alerts for a nitroaniline"
    names = " ".join(a["name"].lower() for a in alerts)
    assert "nitro" in names or "aniline" in names or "amine" in names, names


def test_metformin_clean():
    s = outcome_scores(METFORMIN)
    assert s["hepatotox"]["flagged"] is False, s["hepatotox"]
    assert not reactive_alerts(METFORMIN), reactive_alerts(METFORMIN)


def test_tier_labelled_and_no_probability_keys():
    panel = outcome_panel(TROGLITAZONE)
    assert panel["tier"] == "model-predicted"
    for ep, d in panel["endpoints"].items():
        assert d["tier"] == "model-predicted", (ep, d)
        assert "probability" not in d and "p_toxic" not in d, (ep, d)
