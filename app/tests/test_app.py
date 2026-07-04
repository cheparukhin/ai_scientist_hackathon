"""AppTest render checks (N1 structure images, N3 scope, N4 validation, N6 M2 section).

Run:  .venv/bin/python -m pytest app/tests/test_app.py -q
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from streamlit.testing.v1 import AppTest  # noqa: E402

from core import EXAMPLES  # noqa: E402
from render import mol_png  # noqa: E402
import validation as V  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FINDINGS = os.path.join(REPO, "experiments", "derisk", "expanded", "FINDINGS.md")

APP = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "streamlit_app.py")


def _run(smiles):
    at = AppTest.from_file(APP, default_timeout=60)
    at.session_state["smiles"] = smiles
    at.run()
    return at


# ---------------- N1 ----------------

def test_mol_png_valid_returns_bytes():
    png = mol_png(EXAMPLES["rimonabant"])
    assert isinstance(png, (bytes, bytearray))
    assert len(png) > 0


def test_mol_png_invalid_returns_none():
    assert mol_png("not_a_smiles") is None
    assert mol_png("") is None


def test_app_renders_candidate_structure_image():
    at = _run(EXAMPLES["rimonabant"])
    assert not at.exception
    assert len(at.get("imgs")) >= 1


# ---------------- N3 ----------------

def test_app_has_limits_expander():
    at = _run(EXAMPLES["rimonabant"])
    assert not at.exception
    labels = [e.label for e in at.expander]
    assert any("Limits" in (lbl or "") for lbl in labels), labels


# ---------------- N4 ----------------

def _findings_ranks():
    """Parse the per-drug table in FINDINGS.md -> {drug: (B_rank, ours_rank)} for buried."""
    out = {}
    with open(FINDINGS) as fh:
        for line in fh:
            if not line.strip().startswith("|"):
                continue
            cells = [c.strip().replace("**", "") for c in line.strip().strip("|").split("|")]
            if len(cells) < 6:
                continue
            drug, cat, _culprit, ours, b, _c = cells[0], cells[1], cells[2], cells[3], cells[4], cells[5]
            if cat != "buried":
                continue
            try:
                out[drug] = (int(b), int(ours))
            except ValueError:
                continue
    return out


def test_buried_rank_pairs_match_findings():
    parsed = _findings_ranks()
    assert len(V.BURIED_RANK_PAIRS) == 10
    for drug, default_rank, our_rank in V.BURIED_RANK_PAIRS:
        assert drug in parsed, f"{drug} not found in FINDINGS.md table"
        assert parsed[drug] == (default_rank, our_rank), (drug, parsed[drug], (default_rank, our_rank))


def test_app_validation_panel_present():
    at = _run(EXAMPLES["rimonabant"])
    assert not at.exception
    labels = [e.label for e in at.expander]
    assert any("Track record" in (lbl or "") for lbl in labels), labels
    blob = " ".join(m.value for m in at.markdown)
    assert "11.3" in blob and "3.8" in blob, blob[:400]
    # bar chart: st.bar_chart renders as a vega_lite_chart; its dataset holds one row per
    # (drug, series) -> 10 drugs melted across 2 series = 20 rows, 10 distinct drug indices.
    charts = at.get("vega_lite_chart")
    assert charts, "no bar chart rendered"
    import pyarrow as pa
    reader = pa.ipc.open_stream(charts[0].proto.datasets[0].data.data)
    tbl = reader.read_all()
    idx_col = next(c for c in tbl.column_names if c.startswith("index"))
    n_drugs = len(set(tbl.column(idx_col).to_pylist()))
    assert n_drugs == 10, n_drugs


# ---------------- N6 ----------------

TROGLITAZONE = "Cc1c(C)c2c(c(C)c1O)CCC(C)(COc1ccc(CC3SC(=O)NC3=O)cc1)O2"


def test_app_troglitazone_m2_section_and_m1_headline():
    at = _run(TROGLITAZONE)
    assert not at.exception
    blob = " ".join(m.value for m in at.markdown)
    # lower-confidence liver/metabolism section present and mentions the liver
    assert "Extra checks" in blob, blob[:600]
    assert "liver" in blob.lower(), blob[:600]
    # off-target ranking headline still renders (metric elements present)
    metric_labels = [m.label for m in at.metric]
    assert "We'd run it" in metric_labels, metric_labels


def test_troglitazone_narrative_fallback_mentions_m2_no_probability(monkeypatch):
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    from agent import narrative_report
    from core import score_candidate, build_plan
    from outcome_modules import outcome_panel
    r = score_candidate(TROGLITAZONE)
    plan = build_plan(r)
    m2 = outcome_panel(r["canonical_smiles"])
    txt = narrative_report(r, plan, m2)
    # M2 findings mentioned, but explicitly framed as lower-confidence, and never a probability
    assert "lower-confidence" in txt.lower(), txt
    assert "liver" in txt.lower(), txt
    assert "probability" not in txt.lower(), txt


# ---------------- N7 (stretch) ----------------

def test_mechanism_graph_dot_contents():
    from render import mechanism_graph_dot
    from core import score_candidate, build_plan
    r = score_candidate(EXAMPLES["rimonabant"])
    dot = mechanism_graph_dot(r, build_plan(r))
    assert "digraph" in dot
    assert "CB1" in dot, dot
    assert ("rimonabant" in dot or "taranabant" in dot), dot


def test_app_renders_graphviz_chart():
    at = _run(EXAMPLES["rimonabant"])
    assert not at.exception
    assert len(at.get("graphviz_chart")) >= 1
