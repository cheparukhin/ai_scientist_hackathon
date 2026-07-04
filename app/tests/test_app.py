"""AppTest render checks (N1 structure images, N3 scope, N4 validation, N6 M2 section).

Run:  .venv/bin/python -m pytest app/tests/test_app.py -q
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from streamlit.testing.v1 import AppTest  # noqa: E402

from core import EXAMPLES  # noqa: E402
from render import mol_png  # noqa: E402

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
