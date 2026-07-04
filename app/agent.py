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
    "- If the evidence says the candidate is a known analog of a failed drug, say cheap "
    "similarity already flags it and the engine adds most value for novel chemotypes.\n"
    "- Name the shared mechanism/off-target, the linked failed drug(s), and recommend the "
    "first assay to run and why (severity + resemblance to withdrawn/failed drugs).\n"
    "Write 3-6 sentences, plain prose, no markdown headers, no bullet lists."
)


def _evidence_from(result, plan):
    """Distill result+plan into a compact JSON evidence packet for the LLM."""
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
    return ev


def _fallback(ev):
    """Deterministic templated report (no API key / SDK failure). Fully grounded."""
    if ev.get("status") == "abstain":
        return ("Abstain: " + (ev.get("reason") or "outside applicability domain") + ". "
                "This candidate falls outside the engine's chemical space, so no off-target "
                "assay prioritization is offered. Note the engine is in any case blind to "
                "metabolite-driven and liver toxicity.")
    top = ev["top_liability"]
    parts = []
    if not top["any_liability_flagged"]:
        parts.append(
            f"No strong off-target liability was flagged for this candidate; the highest-priority "
            f"assay is the {top['assay']} at rank {top['our_rank']}, but on weak signal.")
    else:
        drugs = ", ".join(d["name"] for d in ev["linked_failed_drugs"][:3]) or "known failed drugs"
        parts.append(
            f"The candidate most resembles ligands of the target behind {top['organ']}. We move the "
            f"{top['assay']} to rank {top['our_rank']} (default panel rank {top['default_panel_rank']}) "
            f"and tag it '{top['action_tag']}', because the chemotype resembles {drugs}, withdrawn or "
            f"failed for this liability.")
        cite = next((d["citation"] for d in ev["linked_failed_drugs"]), None)
        if cite:
            parts.append(f"Provenance: {cite}")
    if ev["known_analog"]:
        na = ev["nearest_known_failure"]
        parts.append(
            f"Note this candidate is a known analog of {na['name']} (Tanimoto {na['similarity']}), so "
            f"cheap similarity already catches it; the engine adds most value for novel chemotypes.")
    parts.append(
        "Scope: this covers off-target-mediated failures only and is blind to metabolite-driven and "
        "liver toxicity - treat those as weak-coverage, not cleared.")
    return " ".join(parts)


def narrative_report(result, plan=None):
    """Return a grounded med-chemist narrative string. Never raises."""
    if plan is None and result.get("status") == "ok":
        plan = build_plan(result)
    ev = _evidence_from(result, plan)

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
