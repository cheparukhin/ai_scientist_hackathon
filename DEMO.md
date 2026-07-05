# Demo script — Safety-Test Prioritizer

A tight ~4-minute click path. Every number below is what the app actually renders
(deterministic mode, no API key needed). Run it, don't improvise — each beat is chosen.

## Setup (once, before you present)
```bash
.venv/bin/streamlit run app/streamlit_app.py
```
- First load warms the engine (~5s, one-time, shown with a spinner). Do this **before** the room is watching.
- Leave `ANTHROPIC_API_KEY` **unset** — the plain-English summary is then instant and identical every run.
- The page opens already showing the **Rimonabant** result, so there's no blank first screen.

## The one-liner (say this first)
> "Most drug candidates die. The expensive waste is *how long you screen before you find the reason.*
> Paste a molecule and this reorders its safety panel so the test most likely to kill the program runs first —
> by matching it to drugs that actually failed in the clinic for the same off-target reason. It prioritizes
> experiments; it does not claim a molecule is safe."

## Beat 1 — Rimonabant, the headline win  *(already on screen)*
- Point at the verdict: **▶ Run the CB1 counter-screen first** — a standard panel runs it **#15**, we move it to **#1**.
- "This drug was withdrawn for psychiatric/suicidality via CB1. A default panel wouldn't test CB1 until #15."
- Note the yellow **heads-up**: "this is nearly identical to the known failed drug Rimonabant — a simple lookup would already flag it." → *That's the honest setup for the next beat.*

## Beat 2 — Demo mode: the proof point  *(toggle "Demo mode" ON)*
- Green banner: we **hid Rimonabant and its mechanistic cousins** from the database.
- Verdict still lands on **CB1** (z ≈ **+10.8**), and the summary says: *"this is the proof point… recovered from other molecules alone, not from recognising the candidate itself."*
- "So the signal isn't memorized — it's re-derived from *other* CB1 ligands that never failed. This is the novel-chemotype case the tool is built for."

## Beat 3 — Pergolide: consequence beats likelihood  *(click Pergolide, Demo mode OFF)*
- Verdict: **▶ Run the 5-HT2B counter-screen first** — rank **#11 → #2**, the fen-phen valvulopathy mechanism.
- The teaching point: "Pergolide's *strongest* raw resemblance is actually to D2 — its **therapeutic** target. The tool still puts **5-HT2B valvulopathy** first, because it's grounded in drugs that *died* for it. Consequence, weighted over raw similarity."

## Beat 4 — A genuinely novel candidate  *(click Novel candidate)*
- Verdict: a **class-level** match (D2 binder class), and the copy says plainly: *"a class-level structural match, not a link to a specific clinical failure — but still worth screening early."*
- "No known failed drug sits here. The tool tells you that instead of inventing a story — the honesty is the product."

## Beat 5 — The trust beat: it refuses  *(click Cyclosporine)*
- **⛔️ Out of scope.** MW 1203 > 800 — a large macrocycle outside the reference chemical space.
- "It abstains rather than guess. A tool that never says 'I don't know' can't be trusted when it does answer."

## Optional Beat 6 — intellectual honesty  *(Pergolide + Demo mode ON)*
- Only if a judge probes the leave-one-out. The tool reports it **could not** re-derive 5-HT2B once its ergot cousins are hidden — *"the engine's honest limit… it wins on novel chemotypes whose off-target has other known ligands, not on drugs whose only lookalikes are their own withdrawn family."*
- "We show you the failure mode on purpose. Rimonabant recovers cleanly because ChEMBL has many other CB1 ligands; pergolide doesn't, and we say so."

## If asked "how does the number work?"
Open **Show the calculation** under any result: candidate ECFP4 fingerprint → mean top-5 Tanimoto to that
target's ChEMBL binders → z-scored against a fixed 25-drug background of ordinary marketed drugs. **z ≥ 2 flags.**
It's similarity-enrichment, not a probability of harm.

## Track record (the "Track record" expander, top-right)
On buried off-target liabilities, the killer test moves from mean position **#11.3 → #3.8**; top-3 in
**7/10** vs **1/10** for a standard panel; every drug scored with itself hidden. Lower is better.

## Do-not-do list
- Don't demo **Pergolide + Demo mode** as a *win* — it's the honest-limit case (Beat 6 only).
- Don't lead with **Terfenadine/hERG** — it's a validation anchor; standard panels already run hERG first.
- Don't clear the SMILES box to empty — just click an example or paste a new one.
