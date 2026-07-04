#!/usr/bin/env python3
"""
Clean, dependency-light driver for the ProTox-3.0 REST API.
Recovered + reconstructed from the official (Wayback-archived) protox3_api.py.

Endpoints (verified from the official script):
  POST https://tox.charite.de/protox3/src/api_enqueue.php
       data: input_type=smiles|name, input=<SMILES/name>, requested_data=<json list of one space-joined model string>
       -> returns a task_id (text), plus Retry-After header
  POST https://tox.charite.de/protox3/src/api_retrieve.php   data: id=<task_id>  -> 200 when ready, 404 while pending
  Result CSVs (tab-separated):
       https://tox.charite.de/protox3/csv/<task_id>_tox_class.csv     (acute tox: LD50, tox_class)
       https://tox.charite.de/protox3/csv/<task_id>_result.csv        (organ/pathway/MIE/CYP model actives)
       https://tox.charite.de/protox3/csv/<task_id>_tox_targets.csv   (15 Novartis pharmacophore targets, fit 0-3)

Limit: 250 API queries / IP / day. Uses only the Python stdlib (no requests/pandas).

NOTE: This sandbox IP was TLS-dropped by tox.charite.de during the de-risk run
(SSL: UNEXPECTED_EOF). Run this from an un-rate-limited network to collect live data.
"""
import urllib.request, urllib.parse, json, ssl, time, sys, csv, io

ctx = ssl.create_default_context(); ctx.check_hostname = False; ctx.verify_mode = ssl.CERT_NONE
BASE = "https://tox.charite.de/protox3"
ALL_MODELS = ("dili neuro nephro respi cardio carcino immuno mutagen cyto bbb eco clinical nutri "
              "nr_ahr nr_ar nr_ar_lbd nr_aromatase nr_er nr_er_lbd nr_ppar_gamma "
              "sr_are sr_hse sr_mmp sr_p53 sr_atad5 "
              "mie_thr_alpha mie_thr_beta mie_ttr mie_ryr mie_gabar mie_nmdar mie_ampar mie_kar "
              "mie_ache mie_car mie_pxr mie_nadhox mie_vgsc mie_nis "
              "CYP1A2 CYP2C19 CYP2C9 CYP2D6 CYP3A4 CYP2E1")
MODELS = ["acute_tox tox_targets " + ALL_MODELS]   # single space-joined string in a 1-elem list


def post(url, data, timeout=60):
    req = urllib.request.Request(url, data=urllib.parse.urlencode(data).encode())
    return urllib.request.urlopen(req, context=ctx, timeout=timeout)


def enqueue(smiles):
    r = post(BASE + "/src/api_enqueue.php",
             {"input_type": "smiles", "input": smiles, "requested_data": json.dumps(MODELS)})
    tid = r.read().decode().strip()
    wait = int(r.headers.get("Retry-After", 5)) + 1
    return tid, wait


def fetch_csv(task_id, suffix):
    try:
        r = urllib.request.urlopen(BASE + "/csv/%s_%s.csv" % (task_id, suffix), context=ctx, timeout=60)
        return list(csv.reader(io.StringIO(r.read().decode()), delimiter="\t"))
    except Exception:
        return None


def main():
    drugs = json.load(open(sys.argv[1] if len(sys.argv) > 1 else "../drugs.json"))
    out = open("protox_results.tsv", "w")
    out.write("drug\tcategory\tculprit_target\tsection\ttarget\tprediction\tprobability\n")
    for name, d in drugs.items():
        print("Enqueue", name)
        tid, wait = enqueue(d["smiles"]); time.sleep(wait)
        # poll retrieve
        for _ in range(20):
            try:
                post(BASE + "/src/api_retrieve.php", {"id": tid}); break
            except urllib.error.HTTPError as e:
                if e.code == 404: time.sleep(15); continue
                raise
        for suffix, section in [("tox_class", "acute_tox"), ("result", "model"), ("tox_targets", "tox_target")]:
            rows = fetch_csv(tid, suffix)
            if not rows: continue
            for row in rows[1:]:
                cells = row + [""] * 6
                out.write("\t".join([name, d["category"], d["culprit_target"], section,
                                     cells[-3], cells[-2], cells[-1]]) + "\n")
        out.flush(); time.sleep(2)
    out.close(); print("Wrote protox_results.tsv")


if __name__ == "__main__":
    main()
