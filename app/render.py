"""render.py - visualization helpers for the demo UI.

mol_png(smiles) -> PNG bytes | None         (N1: 2D structure depiction)
mechanism_graph_dot(result, plan) -> str    (N7: mechanism-edge network, stretch)

Pure RDKit / Graphviz-DOT; no local graphviz binary needed (the Streamlit frontend
renders DOT). Never raises on bad input - returns None / a minimal DOT instead.
"""
from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")


def mol_png(smiles, size=(320, 240)):
    """2D depiction of a SMILES as PNG bytes. Returns None on unparseable SMILES.

    Uses RDKit rdMolDraw2D.MolDraw2DCairo; falls back to Draw.MolToImage -> PNG bytes.
    """
    if not smiles:
        return None
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    try:
        from rdkit.Chem.Draw import rdMolDraw2D
        d = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
        rdMolDraw2D.PrepareAndDrawMolecule(d, m)
        d.FinishDrawing()
        png = d.GetDrawingText()
        if png:
            return png
    except Exception:
        pass
    # fallback path (PIL-backed)
    try:
        import io
        from rdkit.Chem import Draw
        img = Draw.MolToImage(m, size=size)
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        return buf.getvalue()
    except Exception:
        return None


# ---- organ-system colours for the mechanism graph (N7) ----
_ORGAN_COLOR = {
    "cardiac": "#e57373", "CNS": "#7986cb", "GI": "#81c784",
    "autonomic": "#ffb74d", "other": "#b0bec5",
}


def _organ_system(organ):
    o = (organ or "").lower()
    if any(k in o for k in ("qt", "cardi", "valv", "cav", "nav", "herg", "arrhyth")):
        return "cardiac"
    if any(k in o for k in ("psychiatr", "suicid", "cns", "seroton", "cb1")):
        return "CNS"
    if any(k in o for k in ("colitis", "gi", "ischemic", "retroperiton")):
        return "GI"
    if any(k in o for k in ("autonomic", "muscarin", "cholinerg")):
        return "autonomic"
    return "other"


def mechanism_graph_dot(result, plan):
    """Graphviz DOT: candidate at centre, each flagged target's linked failed drugs as
    nodes, edges labelled by target/mechanism, node colour by organ system."""
    if not plan or plan.get("status") != "ok":
        return "digraph G { label=\"no plan\"; }"

    cand = "candidate"
    lines = [
        "digraph G {",
        "  rankdir=LR;",
        "  node [style=filled, fontname=Helvetica, fontsize=10];",
        f'  "{cand}" [label="candidate", shape=doublecircle, fillcolor="#455a64", fontcolor=white];',
    ]
    seen_drugs = set()
    for row in plan["rows"]:
        if not row.get("flagged") or not row.get("evidence"):
            continue
        tgt = row["assay_name"]
        color = _ORGAN_COLOR[_organ_system(row.get("organ"))]
        tnode = f"t_{row['target_key']}"
        lines.append(f'  "{tnode}" [label="{tgt}", shape=box, fillcolor="{color}"];')
        lines.append(f'  "{cand}" -> "{tnode}" [label="z={row["z"]}"];')
        for e in row["evidence"]:
            dn = e["name"]
            dnode = f"d_{dn}"
            if dnode not in seen_drugs:
                dc = _ORGAN_COLOR[_organ_system(e.get("organ"))]
                lines.append(f'  "{dnode}" [label="{dn}", fillcolor="{dc}"];')
                seen_drugs.add(dnode)
            lines.append(f'  "{tnode}" -> "{dnode}" [label="{e.get("organ", "")}", style=dashed];')
    lines.append("}")
    return "\n".join(lines)
