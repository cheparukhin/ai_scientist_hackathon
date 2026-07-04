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

import json

from core import EXAMPLES, build_plan, score_candidate  # noqa: E402
from agent import narrative_report  # noqa: E402
from render import mol_png  # noqa: E402
import validation as V  # noqa: E402

_REF_FAILURES = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                             "data", "reference_failures.json")))

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

# ---- scope & limits (N3) ----
with st.expander("⚠️ Scope & limits — read before trusting this"):
    st.markdown(
        "- **Covers off-target-mediated failures only.** The engine ranks in-vitro assays by "
        "structural resemblance to ligands of known safety off-targets. It is **blind to "
        "metabolite-driven toxicity and to liver (hepatotox) injury** — those are handled, at a "
        "**lower evidence tier**, by the *model-predicted (M2)* modules below, which are less "
        "validated than the off-target core.\n"
        "- **Narrow addressable population.** Most historical withdrawals are off-panel — on-target "
        "mechanism-based, idiosyncratic hepatotox, or metabolite-driven — which this tool does not "
        "address. The validated value is specifically the *buried off-target* slice.\n"
        "- **No value on hERG / cardiac channels.** Standard panels already front-load hERG at "
        "rank 1, so we add nothing there — hERG is a **validation anchor only**, not a win.\n"
        "- **Advantage is conditional on novelty.** For a candidate that is a close analog of a "
        "known failed drug, cheap similarity already flags it (see the *known-analog* flag); the "
        "engine adds most value for **novel chemotypes**.\n"
        "- **No score here is a probability of harm** — everything is enrichment z-score / "
        "0–100 priority index / rank."
    )

# ---- validation panel (N4) ----
with st.expander("📊 Validation — measured on 20 historical failures (leave-one-out)"):
    st.markdown(
        f"All numbers below are measured, not asserted — cited to "
        f"`{V.FINDINGS_CITATION}` and `{V.SEC4_CITATION}`. Strict leave-one-out: a drug and "
        f"every reference sharing its culprit target are removed before scoring it.\n\n"
        f"**Buried off-target liabilities (n=10) — the validated slice:**\n"
        f"- Mean assays-to-culprit: **{V.MEAN_ASSAYS_TO_CULPRIT_DEFAULT} → "
        f"{V.MEAN_ASSAYS_TO_CULPRIT_OURS}** (default panel → ours)\n"
        f"- Killer assay in top-3: **{V.TOP3_OURS}** vs **{V.TOP3_DEFAULT} (default)**\n"
        f"- **{V.MONEY_WINS} genuine non-obvious wins** (pergolide, cabergoline, methysergide → "
        f"5-HT2B; rimonabant, taranabant → CB1; alosetron → 5-HT3)\n\n"
        f"**Per-target discrimination (AUC, §4):** hERG {V.PER_TARGET_AUC['hERG']}, "
        f"SERT {V.PER_TARGET_AUC['SERT']}, AChE {V.PER_TARGET_AUC['AChE']}, "
        f"MAO-A {V.PER_TARGET_AUC['MAO-A']}.\n\n"
        f"**Scaffold-split:** {V.SCAFFOLD_POOLED} pooled / {V.SCAFFOLD_NOVEL_ISOLATED} "
        f"novel-isolated.  **Ablation:** naive {V.ABLATION_NAIVE} vs R4 {V.ABLATION_R4}."
    )
    st.markdown("**Killer-assay rank — default panel vs ours (10 buried drugs):**")
    rank_df = pd.DataFrame(
        [{"Default panel": d, "Ours": o} for _drug, d, o in V.BURIED_RANK_PAIRS],
        index=[drug for drug, _d, _o in V.BURIED_RANK_PAIRS],
    )
    st.bar_chart(rank_df)
    st.caption("Lower is better (rank 1 = killer assay run first). Fenfluramine / "
               "dexfenfluramine / sibutramine are metabolite-active — the documented boundary "
               "of a parent-structure method.")

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

loo = st.checkbox(
    "Demo mode — leave-one-out (score this known drug as if novel: remove it and its "
    "mechanistic partners from the DB)",
    value=False,
    help="DEMO ONLY. When ON and the input is a known reference-failure drug, the engine "
         "excludes that drug and every reference sharing its culprit target from every "
         "per-target score and from the known-analog check — showing it would still recover "
         "the buried liability on a truly novel chemotype. The live/novel path never removes "
         "anything (default OFF).")

run = st.button("Score candidate", type="primary")

if run or smiles:
    result = score_candidate(smiles, metabolite_smiles=metabolite or None, loo=loo)

    if loo and result.get("loo_matched_ref"):
        st.info(
            f"🧪 **Demo mode active** — scoring **{result['loo_matched_ref']}** as if novel: "
            f"removed it and {len(result['loo_exclude_iks']) - 1} mechanistic partner(s) from "
            f"the reference DB. The recovered liability below comes only from *other* ligands, "
            f"never the drug itself.")
    elif loo:
        st.caption("Demo mode ON, but this input is not a known reference-failure drug — "
                   "leave-one-out only excludes the candidate's own structure.")

    # ---- candidate structure depiction (N1) ----
    cand_png = mol_png(smiles)
    if cand_png:
        st.markdown("### Candidate structure")
        st.image(cand_png, caption="Candidate", width=320)

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

    # ---- side-by-side: candidate vs headline linked failed drug (N1) ----
    head_row = next((r for r in plan["rows"] if r["target_key"] == head["target_key"]), None)
    head_drug = (head_row["evidence"][0]["name"] if head_row and head_row["evidence"] else None)
    if cand_png and head_drug and head_drug in _REF_FAILURES:
        drug_png = mol_png(_REF_FAILURES[head_drug]["smiles"])
        if drug_png:
            st.markdown("#### Candidate vs linked failed drug")
            ic1, ic2 = st.columns(2)
            ic1.image(cand_png, caption="Candidate", width=300)
            ic2.image(drug_png, caption=f"{head_drug} (withdrawn/failed at {head['assay_name']})", width=300)

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
