"""Checks on the panel-2 transmission export.

Deliberately does not import the exporter: the point is to re-derive what can
be re-derived and to cross-check the rest against a record written by a
different script on a different code path.
"""
import base64, csv, json, os, sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "Data", "fig2_studio", "tx")
PREV = os.path.join(ROOT, "tests", "_prevalence_window.csv")
CELLS = ["ll", "ls", "sl", "ss"]

N = {"pass": 0, "fail": 0}


def ck(cond, msg):
    N["pass" if cond else "fail"] += 1
    print(("  ok   " if cond else "  FAIL ") + msg)


def dec(s, n):
    a = np.frombuffer(base64.b64decode(s), np.uint16)
    return a.astype(np.int64) if len(a) == n else None


idx = json.load(open(os.path.join(OUT, "index.json")))
scens = idx["scenarios"]

print("A. index")
ck(len(scens) == 5, f"five scenarios present: {scens}")
ck(idx["cells"] == CELLS, f"cell order is transmitter->recipient: {idx['cells']}")
ck(not idx["failed"], f"no run failed the two-route check ({len(idx['failed'])} failures)")
ck(all(c == 50 for c in idx["counts"].values()), f"50 runs each: {idx['counts']}")

DATA = {s: json.load(open(os.path.join(OUT, f"tx_{s}.json")))["runs"] for s in scens}

print("\nB. every run's own arithmetic closes")
bad_len = bad_sum = bad_tot = bad_neg = 0
for s in scens:
    for r in DATA[s]:
        ser = {c: dec(r[c], r["n_days"]) for c in CELLS}
        if any(v is None for v in ser.values()):
            bad_len += 1
            continue
        if sum(int(v.sum()) for v in ser.values()) != r["n_events"]:
            bad_sum += 1
        if r["n_events"] + r["n_seeded"] != r["n_infected"]:
            bad_tot += 1
        if r["n_long"] + r["n_standard"] != r["n_infected"]:
            bad_neg += 1
ck(bad_len == 0, f"all four series decode to n_days values ({bad_len} bad)")
ck(bad_sum == 0, f"series sum to the recorded event total ({bad_sum} bad)")
ck(bad_tot == 0, f"events + seeded = infected ({bad_tot} bad)")
ck(bad_neg == 0, f"long + standard hosts = infected ({bad_neg} bad)")

print("\nC. host counts agree with the prevalence side-car (separate script)")
if not os.path.exists(PREV):
    ck(False, f"{PREV} missing -- cannot cross-check")
else:
    want = {}
    with open(PREV) as fh:
        for row in csv.DictReader(fh):
            want[(row["scenario"], row["seed"])] = (int(row["n_infected"]),
                                                    int(row["n_long_hosts"]))
    seen = miss = dis = 0
    for s in scens:
        for r in DATA[s]:
            k = (s, r["seed"])
            if k not in want:
                miss += 1
                continue
            seen += 1
            if want[k] != (r["n_infected"], r["n_long"]):
                dis += 1
    ck(seen > 200, f"{seen} runs matched to the side-car")
    ck(miss == 0, f"no run absent from the side-car ({miss} missing)")
    ck(dis == 0, f"infected and long-shedder counts identical on both paths ({dis} differ)")

print("\nD. the day axis")
bad_days = bad_phi = 0
for s in scens:
    for r in DATA[s]:
        if r["n_days"] != int(float(r["end"])) + 1:
            bad_days += 1
        if r["phi"] is not None and not (0 <= float(r["phi"]) <= float(r["end"])):
            bad_phi += 1
ck(bad_days == 0, f"n_days = floor(end)+1 on every run ({bad_days} bad)")
ck(bad_phi == 0, f"burn-in day inside the run on every run that has one ({bad_phi} bad)")
nophi = sum(1 for s in scens for r in DATA[s] if r["phi"] is None)
ck(nophi == 5, f"exactly 5 runs never reach the burn-in, as recorded elsewhere (got {nophi})")

print("\nE. the burn-in cut is a sum range, and it bites")
moved = total = 0
for s in scens:
    for r in DATA[s]:
        if r["phi"] is None:
            continue
        total += 1
        d0 = int(np.floor(float(r["phi"])))
        whole = sum(int(dec(r[c], r["n_days"]).sum()) for c in CELLS)
        after = sum(int(dec(r[c], r["n_days"])[d0:].sum()) for c in CELLS)
        if after < whole:
            moved += 1
        if after > whole:
            moved = -10 ** 9
ck(moved == total, f"post-burn-in total is strictly smaller on all {total} runs (got {moved})")

print("\nF. control has no long shedders, so no long transmitter or recipient")
ctl = DATA["control"]
lv = sum(int(dec(r[c], r["n_days"]).sum()) for r in ctl for c in ("ll", "ls", "sl"))
ck(lv == 0, f"control: zero events involving a long shedder (got {lv})")
ck(all(r["n_long"] == 0 for r in ctl), "control: zero long-shedder hosts")
oth = {s: sum(r["n_long"] for r in DATA[s]) for s in scens if s != "control"}
ck(all(v > 0 for v in oth.values()), f"every other scenario has long-shedder hosts: {oth}")

print(f"\n{N['pass']} passed, {N['fail']} failed")
sys.exit(1 if N["fail"] else 0)
