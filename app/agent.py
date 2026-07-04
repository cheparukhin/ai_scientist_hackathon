"""agent.py - grounded LLM med-chemist narrative.

narrative_report(result, plan=None) -> str

The agent is GROUNDED (BUILD.md sec 2.7): it receives ONLY the retrieved structured
evidence (targets, ranks, linked failed drugs + citations, flags) and is instructed to
invent no drugs, mechanisms, or citations. If ANTHROPIC_API_KEY is unset or the SDK call
fails, it degrades to a deterministic templated report so the app still runs.

Model: claude-opus-4-8 (adaptive thinking). See the claude-api skill for SDK usage.
"""
import json
import os

try:
    from core import build_plan
except ModuleNotFoundError:  # when imported as app.agent
    from app.core import build_plan

MODEL = "claude-opus-4-8"

SYSTEM = (
    "You are a medicinal-chemistry safety advisor writing one short paragraph for a "
    "drug-discovery team. You are given ONLY structured evidence retrieved by a "
    "similarity-based off-target screening engine. Rules you MUST follow:\n"
    "- Use ONLY the drugs, targets, organs, and citations in the evidence. Invent NOTHING "
    "(no drugs, no mechanisms, no citations, no numbers).\n"
    "- Never state a probability of toxicity. The engine outputs enrichment/rank, not P(harm). "
    "Speak in terms of resemblance, priority, and which assay to run first.\n"
    "- The tool sees OFF-TARGET liabilities only. It is blind to metabolite-driven and liver "
    "(hepatotoxicity) toxicity - do not imply the molecule is 'clear'; note this scope limit.\n"
    "- If the evidence says the candidate is a known analog of a failed drug, say a simple "
    "structure lookup already flags it and the tool adds most value on genuinely novel molecules.\n"
    "- Name the shared mechanism/off-target, the linked failed drug(s), and recommend the "
    "first assay to run and why (severity + resemblance to withdrawn/failed drugs).\n"
    "- If model_predicted_modules are present, you MAY add ONE sentence on any liver, "
    "mitochondrial, or reactive-metabolite finding, but ONLY framed as a LOWER-CONFIDENCE "
    "signal that is LESS VALIDATED than the off-target ranking and needs lab confirmation. "
    "Use the phrase 'lower-confidence'. Reactive-metabolite items are structural ALERTS "
    "(hunches), not predictions. Never merge these into the main test-ranking claim and never "
    "state a probability.\n"
    "Write 3-6 sentences, plain prose, no markdown headers, no bullet lists."
)


def _m2_summary(m2):
    """Compact, tier-labelled summary of the M2 model-predicted modules for the packet."""
    if not m2:
        return None
    flagged = {ep: d for ep, d in m2.get("endpoints", {}).items() if d.get("flagged")}
    return {
        "evidence_tier": "lower-confidence (less validated than the off-target ranking; confirm in the lab)",
        "reactive_metabolite_alerts": [
            {"liability": a["name"], "citation": a["citation"], "kind": "structural alert (hypothesis)"}
            for a in m2.get("reactive", [])
        ][:6],
        "flagged_organ_tox": [
            {"endpoint": ep, "enrichment_z": d["z"], "nearest_analog": d["nearest"]["name"],
             "citation": d["nearest"]["citation"]}
            for ep, d in flagged.items()
        ],
        "note": "Similarity-enrichment / structural alerts, never a probability of harm.",
    }


def _evidence_from(result, plan, m2=None):
    """Distill result+plan(+M2) into a compact JSON evidence packet for the LLM."""
    if result.get("status") == "abstain":
        return {
            "status": "abstain",
            "reason": result.get("reason"),
            "message": "The candidate is outside the engine's applicability domain; no off-target "
                       "assay prioritization is offered.",
        }
    head = plan["headline"]
    head_row = next(r for r in plan["rows"] if r["target_key"] == head["target_key"])
    head_drugs = [
        {"name": e["name"], "organ": e["organ"], "tier": e["tier"], "citation": e["citation"]}
        for e in head_row["evidence"]
    ]
    other_drugs = [
        {"name": e["name"], "organ": e["organ"], "citation": e["citation"]}
        for row in plan["rows"] if row["flagged"] and row["target_key"] != head["target_key"]
        for e in row["evidence"]
    ][:4]
    ev = {
        "status": "ok",
        "candidate_smiles": result["canonical_smiles"],
        "top_liability": {
            "assay": head["assay_name"],
            "organ": head["organ"],
            "our_rank": head["our_rank"],
            "default_panel_rank": head["default_rank"],
            "action_tag": head["action"],
            "marginal_value_vs_default_panel": head["marginal_value"],
            "any_liability_flagged": head["any_flagged"],
        },
        "linked_failed_drugs": head_drugs,
        "other_flagged_liabilities": other_drugs,
        "known_analog": result["flags"]["known_analog"],
        "nearest_known_failure": {
            "name": result["nearest_analog"]["name"],
            "similarity": result["nearest_analog"]["sim"],
            "citation": result["nearest_analog"]["citation"],
        },
        "weak_coverage": result["flags"]["weak_coverage"],
        "scope_note": "Off-target-mediated failures only; blind to metabolite and hepatotoxicity.",
    }
    m2s = _m2_summary(m2)
    if m2s:
        ev["model_predicted_modules"] = m2s
    return ev


def _fallback(ev):
    """Deterministic templated report (no API key / SDK failure). Fully grounded."""
    if ev.get("status") == "abstain":
        return ("Out of scope — the tool abstains rather than guess. "
                + (ev.get("reason") or "This molecule is unlike the reference set") + ". "
                "It falls outside the chemical space the tool can judge, so no test ranking is "
                "offered. The tool is in any case blind to metabolite-driven and liver toxicity.")
    top = ev["top_liability"]
    parts = []
    if not top["any_liability_flagged"]:
        parts.append(
            f"No strong off-target liability was flagged for this candidate; the highest-priority "
            f"assay is the {top['assay']} at rank {top['our_rank']}, but on weak signal.")
    else:
        action_plain = {"no-go": "a likely no-go", "counter-screen": "a priority early screen",
                        "monitor": "one to monitor"}.get(top["action_tag"], top["action_tag"])
        linked = ev["linked_failed_drugs"][:3]
        parts.append(
            f"This candidate resembles known binders of the off-target behind {top['organ']}. We move "
            f"the {top['assay']} to #{top['our_rank']} (a standard panel runs it "
            f"#{top['default_panel_rank']}) and mark it {action_plain}.")
        if linked:
            drugs = ", ".join(d["name"].title() for d in linked)
            parts.append(f"Its structure resembles {drugs} — drugs withdrawn or failed for this problem.")
            cite = linked[0].get("citation")
            if cite:
                parts.append(f"Source: {cite}")
        else:
            parts.append(
                "That ranking comes from other molecules linked to this target in the database, not "
                "from recognising the candidate itself.")
    if ev["known_analog"]:
        na = ev["nearest_known_failure"]
        parts.append(
            f"Note this candidate is structurally almost identical to the failed drug {na['name']} "
            f"({int(round((na['similarity'] or 0) * 100))}% similar), so a simple structure lookup "
            f"already flags it; the tool adds most value on genuinely novel molecules.")
    parts.append(
        "Scope: this covers off-target effects only and is blind to metabolite-driven and liver "
        "toxicity - treat those as low-confidence, not cleared.")
    _EP_PLAIN = {"hepatotox": "liver injury", "mito": "mitochondrial toxicity"}
    m2 = ev.get("model_predicted_modules")
    if m2:
        flagged = m2.get("flagged_organ_tox", [])
        alerts = m2.get("reactive_metabolite_alerts", [])
        if flagged:
            eps = ", ".join(f"{_EP_PLAIN.get(d['endpoint'], d['endpoint'])} "
                            f"(closest: {d['nearest_analog'].title()})" for d in flagged)
            parts.append(
                f"Separately, a lower-confidence check — less validated than the ranking above, so "
                f"confirm it in the lab — flags possible {eps}.")
        if alerts:
            liabs = ", ".join(a["liability"] for a in alerts[:3])
            parts.append(
                f"It also raises reactive-metabolite structural alerts (hunches to confirm, not "
                f"predictions): {liabs}.")
        if not flagged and not alerts:
            parts.append(
                "The lower-confidence liver/metabolism check flagged nothing — but that is still "
                "not a clean bill of health.")
    return " ".join(parts)


def narrative_report(result, plan=None, m2=None):
    """Return a grounded med-chemist narrative string. Never raises.

    m2: optional outcome_modules.outcome_panel(...) dict - mentioned ONLY as a lower-tier,
    model-predicted signal needing confirmation (never merged into the M1 claim).
    """
    if plan is None and result.get("status") == "ok":
        plan = build_plan(result)
    ev = _evidence_from(result, plan, m2)

    api_key = os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        return _fallback(ev)

    try:
        import anthropic

        client = anthropic.Anthropic(api_key=api_key)
        user = (
            "Here is the retrieved structured evidence (JSON). Write the grounded paragraph "
            "described in your instructions, using only what is present here:\n\n"
            + json.dumps(ev, indent=1)
        )
        resp = client.messages.create(
            model=MODEL,
            max_tokens=1024,
            thinking={"type": "adaptive"},
            system=SYSTEM,
            messages=[{"role": "user", "content": user}],
        )
        text = "".join(b.text for b in resp.content if b.type == "text").strip()
        return text or _fallback(ev)
    except Exception as e:  # graceful degradation - the app must still run
        return _fallback(ev) + f"\n\n_(LLM narrative unavailable: {type(e).__name__}; showed templated report.)_"
