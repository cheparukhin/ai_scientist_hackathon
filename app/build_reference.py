"""Build app/data/reference_failures.json from experiments/derisk/drugs.json plus
curated provenance/citations (toxicity-assay-recommender.md sec 9).

Every reference drug MUST carry a citation (honesty guardrail sec 5). Tier T1 =
clinical failure / withdrawn / restricted (these 20 are all program-ending
liabilities used as the grounding set).

Run:  .venv/bin/python app/build_reference.py
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DRUGS = json.load(open(os.path.join(HERE, "..", "experiments", "derisk", "drugs.json")))

# Curated provenance per drug (grounded; no invented citations). Sources: toxicity-
# assay-recommender.md sec 9 "Verified demo mechanisms - citations" + regulatory record.
CITATIONS = {
    "pergolide":       "Zanettini et al., NEJM 2007 (PMID 17202454); Setola et al., Mol Pharmacol 2003 (PMID 12761331). Withdrawn (US 2007) - valvular heart disease.",
    "cabergoline":     "Zanettini et al., NEJM 2007 (PMID 17202454). Dose-restricted - cardiac valvulopathy.",
    "fenfluramine":    "Connolly et al., NEJM 1997 (PMID 9271479). Withdrawn 1997 - valvulopathy (metabolite-active).",
    "dexfenfluramine": "Connolly et al., NEJM 1997 (PMID 9271479). Withdrawn 1997 - valvulopathy (metabolite-active).",
    "methysergide":    "Setola et al., Mol Pharmacol 2003 (PMID 12761331). Retroperitoneal/valve fibrosis.",
    "rimonabant":      "Christensen et al., Lancet 2007 (PMID 18029027). Withdrawn (EU 2008) - psychiatric / suicidality.",
    "taranabant":      "Aronne et al., 2010. Development discontinued - psychiatric adverse events (CB1).",
    "alosetron":       "FDA 2000 withdrawal / restriction - ischemic colitis (5-HT3).",
    "sibutramine":     "SCOUT trial; FDA & EMA 2010 withdrawal - cardiovascular events (NET).",
    "mibefradil":      "Withdrawn 1998 - fatal CYP3A4 drug-drug interactions + cardiac (Cav1.2 / T-type).",
    "terfenadine":     "Saxena et al., PMC2845965. Withdrawn 1997 - QT / torsades (hERG).",
    "cisapride":       "Saxena et al., PMC2845965. Withdrawn 2000 - QT / torsades (hERG).",
    "astemizole":      "Withdrawn 1999 - QT / torsades (hERG).",
    "grepafloxacin":   "Withdrawn 1999 - QT / torsades (hERG).",
    "sertindole":      "Suspended 1998 (EU) - QT / arrhythmia (hERG).",
    "thioridazine":    "Restricted / withdrawn - QT prolongation (hERG).",
    "terodiline":      "Withdrawn 1991 - QT / torsades (hERG).",
    "droperidol":      "FDA black-box 2001 - QT prolongation (hERG).",
    "sparfloxacin":    "Withdrawn - QT prolongation + phototoxicity (hERG).",
    "mesoridazine":    "Withdrawn 2004 - QT prolongation (hERG).",
}

# culprit_target string -> panel target key (matches panel.json keys)
CULPRIT_TO_KEY = {
    "5-HT2B": "5HT2B_CHEMBL1833",
    "CB1":    "CB1_CHEMBL218",
    "5-HT3":  "5HT3_CHEMBL1899",
    "NET":    "NET_CHEMBL222",
    "Cav1.2": "Cav12_CHEMBL1940",
    "hERG":   "hERG_CHEMBL240",
}

out = {}
missing = []
for name, d in DRUGS.items():
    if name not in CITATIONS:
        missing.append(name)
        continue
    out[name] = {
        "name": name,
        "smiles": d["smiles"],
        "culprit_target": d["culprit_target"],
        "culprit_target_key": CULPRIT_TO_KEY[d["culprit_target"]],
        "culprit_chembl_id": d.get("culprit_chembl_id", ""),
        "organ": d["organ"],
        "category": d["category"],
        "tier": "T1",
        "citation": CITATIONS[name],
    }

if missing:
    raise SystemExit(f"Missing citations for: {missing} - refusing to write (every ref drug needs provenance).")

path = os.path.join(HERE, "data", "reference_failures.json")
json.dump(out, open(path, "w"), indent=1)
print(f"Wrote {path}: {len(out)} reference failed drugs, all with citations.")
for k, v in out.items():
    print(f"  {k:16s} {v['culprit_target']:8s} {v['tier']}  {v['citation'][:60]}...")
