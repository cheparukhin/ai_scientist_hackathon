"""streamlit_app.py - the demo UI.

Paste a small-molecule SMILES (or click an example) -> reordered in-vitro safety-assay
plan that moves the assay most likely to catch a program-ending off-target liability to
the front, with the assays-to-culprit rank delta, an evidence trail of linked failed
drugs + citations, and a grounded LLM med-chemist narrative.

Run:  .venv/bin/streamlit run app/streamlit_app.py
"""
import os
import sys

import pandas as pd
import streamlit as st

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from core import EXAMPLES, build_plan, score_candidate  # noqa: E402
from agent import narrative_report  # noqa: E402

st.set_page_config(page_title="Toxicity Assay Recommender", page_icon="🧪", layout="wide")


@st.cache_resource(show_spinner="Loading 18-target panel + reference DB (one-time)...")
def _warm_engine():
    from core import get_engine
    get_engine()
    return True


_warm_engine()

# ---- example molecules (label -> smiles); note which is which for the honest framing ----
EXAMPLE_BUTTONS = [
    ("Rimonabant (CB1 - buried)", EXAMPLES["rimonabant"]),
    ("Pergolide (5-HT2B - buried)", EXAMPLES["pergolide"]),
    ("Terfenadine (hERG - anchor)", EXAMPLES["terfenadine"]),
    ("Novel aryl-piperazine", EXAMPLES["novel_arylpiperazine"]),
    ("Cyclosporine (abstain)", EXAMPLES["cyclosporine"]),
]

ACTION_EMOJI = {"no-go": "🛑 no-go", "counter-screen": "🔬 counter-screen", "monitor": "👁 monitor"}

st.title("🧪 Toxicity Risk → Assay Recommender")
st.caption(
    "Reorders an in-vitro safety panel to run the **killer** off-target counter-screen first. "
    "Scores are enrichment / rank vs a background set - **never a probability of harm**. "
    "Covers *off-target-mediated* failures only; **blind to metabolite-driven and liver "
    "(hepatotox) toxicity** - those are weak-coverage, not 'clear'."
)

# ---- input ----
if "smiles" not in st.session_state:
    st.session_state.smiles = EXAMPLES["rimonabant"]

st.markdown("**Examples:**")
cols = st.columns(len(EXAMPLE_BUTTONS))
for col, (label, smi) in zip(cols, EXAMPLE_BUTTONS):
    if col.button(label, width="stretch"):
        st.session_state.smiles = smi

smiles = st.text_input("Candidate SMILES", key="smiles")
metabolite = st.text_input(
    "Known active metabolite SMILES (optional)", value="",
    help="If given, parent AND metabolite are scored and aggregated by MAX per target.")

run = st.button("Score candidate", type="primary")

if run or smiles:
    result = score_candidate(smiles, metabolite_smiles=metabolite or None)

    # ---------------- ABSTAIN ----------------
    if result.get("status") == "abstain":
        st.error(
            f"**Abstain — outside applicability domain (gate {result.get('gate')}).**  \n"
            f"{result.get('reason')}"
        )
        st.info(
            "The descriptor-box scope gate rejects molecules unlike the reference set (large "
            "peptides/macrocycles, highly halogenated, sugars/polyols, metals, or too small). "
            "No off-target assay prioritization is offered for out-of-domain candidates."
        )
        with st.spinner("Generating narrative..."):
            st.markdown("#### Med-chemist narrative")
            st.write(narrative_report(result))
        st.stop()

    # ---------------- OK ----------------
    plan = build_plan(result)
    head = plan["headline"]

    # banners
    if not head["any_flagged"]:
        st.warning(
            "No strong off-target liability flagged (best signal is weak). The plan below is "
            "low-confidence; this is **weak-coverage**, not a clean bill of health."
        )
    if result["flags"]["known_analog"]:
        na = result["nearest_analog"]
        st.warning(
            f"**Known-analog flag:** highly similar to failed drug **{na['name']}** "
            f"(Tanimoto {na['sim']}). Cheap similarity already catches this — the engine adds "
            f"most value for *novel* chemotypes."
        )
    if head["target_key"] == "hERG_CHEMBL240":
        st.info("Top hit is hERG — **standard panels already front-load this**, so the marginal "
                "value of reordering is low here. The engine's edge is on buried off-targets.")

    # headline metric: assays-to-culprit
    st.markdown("### Assays-to-culprit")
    c1, c2, c3 = st.columns(3)
    c1.metric("Killer assay", head["assay_name"])
    c2.metric("Our rank", f"#{head['our_rank']}")
    c3.metric("Default-panel rank", f"#{head['default_rank']}", delta=f"{head['delta']:+d} moved earlier",
              delta_color="normal")
    st.caption(f"Organ / phenotype: {head['organ']}  ·  recommended action: **{head['action']}**  ·  "
               f"marginal value vs default panel: **{head['marginal_value']}**")

    # ---- reordered plan table ----
    st.markdown("### Reordered assay plan")
    rows = plan["rows"]
    df = pd.DataFrame([{
        "Our rank": r["our_rank"],
        "Assay": r["assay_name"],
        "Action": ACTION_EMOJI.get(r["action"], r["action"]),
        "z-score": r["z"],
        "Priority (0-100)": r["priority"],
        "Default rank": r["default_rank"],
        "Δ (earlier)": r["delta"],
        "Organ": r["organ"],
        "Severity": r["severity"],
        "Flagged": "●" if r["flagged"] else "",
    } for r in rows])
    st.dataframe(df, hide_index=True, width="stretch",
                 column_config={"Δ (earlier)": st.column_config.NumberColumn(format="%+d")})

    # ---- evidence trail ----
    st.markdown("### Evidence trail — linked failed drugs")
    flagged = [r for r in rows if r["flagged"] and r["evidence"]]
    if not flagged:
        st.write("No flagged target carries known-failure grounding for this candidate.")
    for r in flagged:
        st.markdown(f"**{r['assay_name']}** · {r['organ']} · _{r['action']}_")
        ev_df = pd.DataFrame([{
            "Failed drug": e["name"], "Organ / phenotype": e["organ"],
            "Tier": e["tier"], "Provenance / citation": e["citation"],
        } for e in r["evidence"]])
        st.dataframe(ev_df, hide_index=True, width="stretch")

    # ---- LLM narrative ----
    st.markdown("### Med-chemist narrative")
    key_set = bool(os.environ.get("ANTHROPIC_API_KEY"))
    with st.spinner("Generating grounded narrative..."):
        st.write(narrative_report(result, plan))
    if not key_set:
        st.caption("_ANTHROPIC_API_KEY not set — showing the deterministic templated report "
                   "(fully grounded). Set the key for the LLM-generated version._")

    # metabolite note
    if result.get("metabolite"):
        st.caption(f"Metabolite handling: {result['metabolite'].get('aggregation', result['metabolite'])}")
