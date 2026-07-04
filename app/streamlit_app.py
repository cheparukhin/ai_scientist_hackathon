"""streamlit_app.py - the demo UI.

Plain-language front end: paste a drug candidate, get back the ONE safety test worth
running first (and the full reordered plan), with the real failed drugs it resembles as
evidence, a lower-confidence liver/metabolism check, and a plain-English summary.

Internal terminology (M1/M2/R4, "assays-to-culprit", z-scores, "known-analog") is kept OUT
of the interface - it lives in the code and docs, not in front of a judge.

Run:  .venv/bin/streamlit run app/streamlit_app.py
"""
import json
import os
import sys

import pandas as pd
import streamlit as st

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from core import EXAMPLES, build_plan, score_candidate  # noqa: E402
from agent import narrative_report  # noqa: E402
from render import mol_png, mechanism_graph_dot  # noqa: E402
from outcome_modules import outcome_panel  # noqa: E402
import validation as V  # noqa: E402

_REF_FAILURES = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                             "data", "reference_failures.json")))

st.set_page_config(page_title="Safety-Test Prioritizer", page_icon="🧪", layout="wide")


@st.cache_resource(show_spinner="Loading the safety panel + failed-drug database (one-time)…")
def _warm_engine():
    from core import get_engine
    get_engine()
    return True


_warm_engine()

# ---- example molecules: plain-English organ hint, no jargon ----
EXAMPLE_BUTTONS = [
    ("Rimonabant", EXAMPLES["rimonabant"], "psychiatric"),
    ("Pergolide", EXAMPLES["pergolide"], "heart valve"),
    ("Terfenadine", EXAMPLES["terfenadine"], "heart rhythm"),
    ("Novel candidate", EXAMPLES["novel_arylpiperazine"], "unknown"),
    ("Cyclosporine", EXAMPLES["cyclosporine"], "out of scope"),
]

# plain-language action labels (internal: no-go / counter-screen / monitor)
ACTION_PLAIN = {
    "no-go":          "🛑 Likely no-go",
    "counter-screen": "🔬 Screen early",
    "monitor":        "👁 Monitor",
}
ENDPOINT_PLAIN = {"hepatotox": "Liver injury", "mito": "Mitochondrial"}

# ======================================================================================
#  HEADER
# ======================================================================================
st.title("🧪 Which safety test should you run first?")
st.markdown(
    "Paste a drug candidate. This tool **reorders its safety tests** so the one most likely "
    "to end the program runs first — by matching it to real drugs that were **withdrawn or "
    "failed in the clinic** for the same off-target reason. It tells you *what to run first*, "
    "not whether a molecule is safe."
)

# collapsed trust/limits — visible but out of the way
lc, rc = st.columns(2)
with lc:
    with st.expander("Limits — what this tool can't see"):
        st.markdown(
            "- **It only sees off-target effects** — a candidate binding a protein it shouldn't. "
            "It is **blind to liver injury and metabolite-driven toxicity**; those get a "
            "separate, **lower-confidence** check further down, not the main ranking.\n"
            "- **No value on the heart-rhythm (hERG) test** — every standard panel already runs "
            "that first, so we add nothing there. Our edge is the *buried* off-targets panels "
            "run late.\n"
            "- **Best on genuinely new molecules.** If a candidate is nearly identical to a known "
            "bad drug, a simple lookup already catches it — we earn our keep on novel structures.\n"
            "- **Nothing here is a probability of harm.** Everything is a *match-strength* signal "
            "and a *ranking*."
        )
with rc:
    with st.expander("Track record — measured on 20 real drug failures"):
        st.markdown(
            f"On the buried off-target liabilities standard panels run late, the killer test "
            f"moves from an average position of **#{V.MEAN_ASSAYS_TO_CULPRIT_DEFAULT} → "
            f"#{V.MEAN_ASSAYS_TO_CULPRIT_OURS}** — you reach the go/no-go decision in about "
            f"**{V.MEAN_ASSAYS_TO_CULPRIT_OURS} experiments instead of "
            f"{V.MEAN_ASSAYS_TO_CULPRIT_DEFAULT}**. The killer test lands in the top-3 for "
            f"**{V.TOP3_OURS}** candidates vs **{V.TOP3_DEFAULT}** for a standard panel, with "
            f"**{V.MONEY_WINS} wins a standard panel misses entirely**.\n\n"
            f"_Every drug was scored with itself hidden from the database, so nothing is "
            f"memorised; cited to `{V.FINDINGS_CITATION}`. Bars show where the killer test "
            f"lands — lower is better._"
        )
        rank_df = pd.DataFrame(
            [{"Standard panel": d, "Ours": o} for _drug, d, o in V.BURIED_RANK_PAIRS],
            index=[drug for drug, _d, _o in V.BURIED_RANK_PAIRS],
        )
        st.bar_chart(rank_df)

st.divider()

# ======================================================================================
#  INPUT
# ======================================================================================
if "smiles" not in st.session_state:
    st.session_state.smiles = EXAMPLES["rimonabant"]

st.markdown("##### Try an example")
cols = st.columns(len(EXAMPLE_BUTTONS))
for col, (label, smi, hint) in zip(cols, EXAMPLE_BUTTONS):
    if col.button(f"{label}", width="stretch", help=f"failed for: {hint}"):
        st.session_state.smiles = smi
    col.caption(f"<div style='text-align:center'>{hint}</div>", unsafe_allow_html=True)

smiles = st.text_input("…or paste a SMILES", key="smiles")

loo = st.toggle(
    "Demo mode — score a known drug as if it were brand-new",
    value=False,
    help="When ON and you enter a known failed drug, the tool hides that drug (and its "
         "mechanistic cousins) from its database, then shows it would STILL flag the right "
         "test — proving the result doesn't come from recognising the drug itself. Leave OFF "
         "for a real, novel candidate.")

with st.expander("Advanced — add a known active metabolite"):
    metabolite = st.text_input(
        "Active metabolite SMILES (optional)", value="",
        help="If the real culprit is a metabolite, paste it here; the parent and metabolite "
             "are scored together (whichever looks worse at each target wins).")

run = st.button("Prioritize tests", type="primary")

# ======================================================================================
#  RESULT
# ======================================================================================
if run or smiles:
    result = score_candidate(smiles, metabolite_smiles=metabolite or None, loo=loo)

    st.divider()

    # candidate structure (always, if parseable)
    cand_png = mol_png(smiles)
    top_l, top_r = st.columns([1, 2])
    with top_l:
        if cand_png:
            st.image(cand_png, caption="Your candidate", width=260)
        else:
            st.caption("(structure could not be drawn)")

    # ---------------- OUT OF SCOPE (abstain) ----------------
    if result.get("status") == "abstain":
        with top_r:
            st.subheader("⛔️ Out of scope — no ranking offered")
            st.markdown(
                f"This molecule is unlike anything in our reference set, so we **decline to "
                f"guess** rather than give you a misleading answer.\n\n**Why:** {result.get('reason')}"
            )
            st.caption("The tool only judges ordinary organic small molecules. Very large "
                       "molecules, peptides/macrocycles, sugars, and metal-containing compounds "
                       "are deliberately refused.")
        with st.spinner("Writing summary…"):
            st.markdown("##### Plain-English summary")
            st.write(narrative_report(result))
        st.stop()

    # ---------------- OK ----------------
    plan = build_plan(result)
    head = plan["headline"]
    head_row = next((r for r in plan["rows"] if r["target_key"] == head["target_key"]), None)
    flagged_head = bool(head["any_flagged"])
    is_herg = head["target_key"] == "hERG_CHEMBL240"

    # demo-mode confirmation
    if loo and result.get("loo_matched_ref"):
        st.success(
            f"**Demo mode:** scoring **{result['loo_matched_ref'].title()}** as if it were new — "
            f"we hid it and {len(result['loo_exclude_iks']) - 1} mechanistic cousin(s) from the "
            f"database. Everything below is recovered from *other* molecules, never the drug itself."
        )
    elif loo:
        st.caption("Demo mode is on, but this isn't a known failed drug — nothing extra was hidden.")

    # ---------------- VERDICT CARD ----------------
    with top_r:
        organ = head["organ"]
        if flagged_head and not is_herg:
            st.subheader(f"▶ Run the {head['assay_name']} first")
            st.markdown(
                f"This candidate looks like drugs that failed for **{organ}**. A standard panel "
                f"runs that test **#{head['default_rank']}** — we move it to **#{head['our_rank']}**."
            )
        elif flagged_head and is_herg:
            st.subheader(f"▶ {head['assay_name']} (already standard)")
            st.markdown(
                "The strongest signal is the **heart-rhythm (hERG) test** — but every panel "
                "already runs that first, so there's little to reorder here. This tool's real "
                "value is the *buried* off-targets, not hERG."
            )
        else:
            st.subheader("No strong off-target red flag")
            st.markdown(
                f"Nothing scores as a clear off-target risk (the strongest is the "
                f"**{head['assay_name']}**, on a weak signal). Treat this as **low confidence, "
                f"not a clean bill of health** — and check the liver/metabolism section below."
            )

        m1, m2c, m3 = st.columns(3)
        m1.metric("Top-priority test", head["assay_name"])
        m2c.metric("Our order", f"#{head['our_rank']}")
        delta = head["default_rank"] - head["our_rank"]
        m3.metric("Standard panel order", f"#{head['default_rank']}",
                  delta=(f"+{delta} tests earlier" if delta > 0 else "already up front"),
                  delta_color="normal")
        if flagged_head and not is_herg and delta > 0:
            st.markdown(
                f"➜ **≈{delta} fewer experiments** before you'd reach the go/no-go decision.  "
                f"Recommended call: **{ACTION_PLAIN.get(head['action'], head['action'])}**."
            )

    # known-drug heads-up (only when we did NOT hide it) — teaches why demo mode exists
    if result["flags"]["known_analog"] and not loo:
        na = result["nearest_analog"]
        st.warning(
            f"**Heads-up:** this candidate is nearly identical to the known failed drug "
            f"**{na['name'].title()}**, so a simple structure lookup would already flag it. "
            f"To see the tool work on a *novel* molecule, turn on **Demo mode** above — it hides "
            f"{na['name'].title()} and re-derives the answer from scratch."
        )

    # ---------------- WHY: the failed drugs it resembles ----------------
    st.markdown("### Why — the real failed drugs it resembles")
    head_evidence = head_row["evidence"] if (head_row and flagged_head) else []
    if head_evidence:
        st.markdown(
            f"Its structure matches known binders of the off-target behind **{head['organ']}**. "
            f"Drugs that hit this target and **failed or were withdrawn**:"
        )
        ev_df = pd.DataFrame([{
            "Failed drug": e["name"].title(),
            "What happened": e["organ"],
            "Source": e["citation"],
        } for e in head_evidence])
        st.dataframe(ev_df, hide_index=True, width="stretch")

        # side-by-side structures (candidate vs the closest failed drug)
        head_drug = head_evidence[0]["name"]
        if cand_png and head_drug in _REF_FAILURES:
            drug_png = mol_png(_REF_FAILURES[head_drug]["smiles"])
            if drug_png:
                ic1, ic2 = st.columns(2)
                ic1.image(cand_png, caption="Your candidate", width=280)
                ic2.image(drug_png, caption=f"{head_drug.title()} — failed for {head_evidence[0]['organ']}",
                          width=280)
    elif flagged_head and loo and result.get("loo_matched_ref"):
        st.success(
            f"**This is the proof point.** We hid every known failed drug for this target, yet "
            f"the tool still ranked the **{head['assay_name']}** first — recovered purely from "
            f"other molecules in the database, not from recognising "
            f"{result['loo_matched_ref'].title()} itself."
        )
    else:
        st.caption("No known failed drug is closely linked to the top-ranked test for this candidate.")

    # other flagged tests with evidence -> tucked away
    others = [r for r in plan["rows"]
              if r["flagged"] and r["evidence"] and r["target_key"] != head["target_key"]]
    if others:
        with st.expander(f"See {len(others)} more flagged test(s) and their linked failed drugs"):
            for r in others:
                st.markdown(f"**{r['assay_name']}** · {r['organ']}")
                st.dataframe(pd.DataFrame([{
                    "Failed drug": e["name"].title(), "What happened": e["organ"],
                    "Source": e["citation"],
                } for e in r["evidence"]]), hide_index=True, width="stretch")

    # ---------------- FULL REORDERED PLAN ----------------
    st.markdown("### The full reordered test plan")
    st.caption("Every test in the panel, reordered. **#1 = run first.** "
               "*Match strength* is how strongly the candidate resembles known bad actors at that "
               "target (higher = stronger; not a probability).")
    rows = plan["rows"]
    df = pd.DataFrame([{
        "#": r["our_rank"],
        "Safety test": r["assay_name"],
        "What to do": ACTION_PLAIN.get(r["action"], r["action"]),
        "Match strength": r["z"],
        "Standard order": r["default_rank"],
        "Moved earlier": r["delta"],
        "Organ / effect": r["organ"],
    } for r in rows])
    st.dataframe(
        df, hide_index=True, width="stretch",
        column_config={
            "Match strength": st.column_config.NumberColumn(
                help="Std deviations above a typical drug. Higher = stronger match to known "
                     "bad actors. Not a probability of harm.", format="%.1f"),
            "Moved earlier": st.column_config.NumberColumn(
                help="Positions moved up vs a standard panel.", format="%+d"),
        })

    # mechanism map — optional, tucked in an expander so it never clutters the main flow
    if any(r["flagged"] and r["evidence"] for r in rows):
        with st.expander("See the mechanism map (candidate → off-target → failed drugs)"):
            try:
                st.graphviz_chart(mechanism_graph_dot(result, plan))
            except Exception as e:
                st.caption(f"(map unavailable: {type(e).__name__})")

    # ---------------- LOWER-CONFIDENCE: liver / metabolism ----------------
    st.divider()
    m2 = outcome_panel(result["canonical_smiles"])
    st.markdown("### Extra checks: liver, mitochondria & reactive metabolites")
    st.caption(
        "A separate, **lower-confidence** look at the things the main off-target ranking can't "
        "see. This is *not* validated the way the ranking above is — read it as a prompt for "
        "more lab work, never as a verdict, and never a probability of harm."
    )

    any_liver_flag = any(d["flagged"] for d in m2["endpoints"].values())
    lcol, rcol = st.columns(2)

    with lcol:
        st.markdown("**Organ-toxicity look-alikes**")
        ep_df = pd.DataFrame([{
            "Concern": ENDPOINT_PLAIN.get(ep, ep),
            "Signal": d["z"],
            "Flag": "⚠️ elevated" if d["flagged"] else "—",
            "Most similar known-toxic drug": f"{d['nearest']['name'].title()} "
                                             f"({int(d['nearest']['sim'] * 100)}% similar)",
        } for ep, d in m2["endpoints"].items()])
        st.dataframe(ep_df, hide_index=True, width="stretch",
                     column_config={"Signal": st.column_config.NumberColumn(format="%.1f")})

    with rcol:
        st.markdown("**Reactive-metabolite alerts** — structural hunches to confirm in the lab")
        if m2["reactive"]:
            for a in m2["reactive"][:6]:
                st.markdown(f"- ⚠️ **{a['name']}**")
            if len(m2["reactive"]) > 6:
                st.caption(f"…and {len(m2['reactive']) - 6} more.")
        else:
            st.markdown("- ✅ No reactive-metabolite alert matched.")

    if not any_liver_flag and not m2["reactive"]:
        st.info("Nothing flagged here — but that's **not** a clean bill of health; this check is "
                "the least validated part of the tool.")

    # ---------------- PLAIN-ENGLISH SUMMARY ----------------
    st.markdown("### Plain-English summary")
    key_set = bool(os.environ.get("ANTHROPIC_API_KEY"))
    with st.spinner("Writing summary…"):
        st.write(narrative_report(result, plan, m2))
    if not key_set:
        st.caption("_Auto-generated from the evidence above (deterministic). "
                   "Set ANTHROPIC_API_KEY for the LLM-written version._")

    if result.get("metabolite") and result["metabolite"].get("aggregation"):
        st.caption("A metabolite was included — parent and metabolite were scored together.")
