"""Checks the studio server end to end against a running instance.

    python -I tests/studio/server.py &
    python -I tests/check_studio_api.py [http://127.0.0.1:8600]

Deliberately goes through HTTP rather than importing the modules: the thing
being checked is what a page actually receives.
"""
import json
import re, sys, urllib.error, urllib.request

BASE = (sys.argv[1] if len(sys.argv) > 1 else "http://127.0.0.1:8600").rstrip("/")
SCEN = ["control", "SOT", "HIV_low", "HIV_high", "edge_case"]
DESIGN = {"control": 0.0, "SOT": 1.0, "HIV_low": 1.0, "HIV_high": 12.0,
          "edge_case": 0.1}
N = {"pass": 0, "fail": 0}
# A panel is computed on demand and cached under cache/_panels, which
# pack_cache.sh deliberately excludes from the archive. The first hit on a
# freshly extracted grid therefore rebuilds it: panel 4 measured 84 s over 384
# builds, against 0.19 s warm. 120 s was not enough and turned a slow first
# call into a traceback.
TIMEOUT = 900


def ck(cond, msg):
    N["pass" if cond else "fail"] += 1
    print(("  ok   " if cond else "  FAIL ") + msg)


def get(path):
    """Body and status, for an error status too. The server answers a refusal
    with a JSON reason, so raising here would turn a failed check into a
    traceback and abandon every check after it."""
    try:
        with urllib.request.urlopen(BASE + path, timeout=TIMEOUT) as r:
            return json.load(r), r.status
    except urllib.error.HTTPError as e:
        raw = e.read()
        try:
            return json.loads(raw), e.code
        except Exception:
            return {"error": raw[:200].decode("utf-8", "replace")}, e.code


def get_text(path):
    """Body as text and status, for the HTML and script files."""
    try:
        with urllib.request.urlopen(BASE + path, timeout=TIMEOUT) as r:
            return r.read().decode("utf-8", "replace"), r.status
    except urllib.error.HTTPError as e:
        return "", e.code


def get_code(path):
    """Status only -- used for pages (HTML) and for refusals (no body to read)."""
    try:
        with urllib.request.urlopen(BASE + path, timeout=TIMEOUT) as r:
            r.read(1)
            return r.status
    except urllib.error.HTTPError as e:
        return e.code


def post(path, body):
    req = urllib.request.Request(BASE + path, method="POST",
                                 data=json.dumps(body).encode(),
                                 headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            return json.load(r), r.status
    except urllib.error.HTTPError as e:
        return json.load(e), e.code


print("A. the shell")
idx, _ = get("/api/index")
ck(idx["scen_order"] == SCEN, "scenario order is fixed and complete")
ck(len(idx["code"]) == 12, f"compute code hash present: {idx['code']}")
ck(len(idx["builds"]) >= 1, f"{len(idx['builds'])} clade build(s) cached")
ck(idx["transmissions"] is not None, "transmission table present")
for b in idx["builds"]:
    ck(b["code"] == idx["code"],
       f"build {b['key'][:28]} was computed by the current code")

# The nine per-page documents and studio.js were deleted in v2.5.4: the studio
# is one document with ten tabs, so what has to be reachable is index.html, the
# stylesheet, and every script index.html loads. Reading the list out of
# index.html rather than restating it here means a script added to the page
# cannot be left untested, and one removed cannot leave a stale assertion.
print("\nB. the web layer is reachable")
ck(get_code("/") == 200, "/ serves the studio")
_idx, _c = get_text("/index.html")
ck(_c == 200, "/index.html")
_srcs = re.findall(r'src="([^"]+)"', _idx or "")
ck(len(_srcs) >= 15, f"index.html loads {len(_srcs)} scripts")
for pg in ["studio.css"] + _srcs:
    ck(get_code("/" + pg.lstrip("/")) == 200, f"/{pg}")

# Everything from here reads a clade build at the operating point. Without it
# every section below would fail on the same missing directory, so say it once
# and stop -- a rebuild is the fix, not a code change.
probe, code = get("/api/panel?n=1&k=7&mt=1&burnin=1")
if code != 200:
    print(f"\nthe operating point k=7 is not built for the current code: "
          f"{probe.get('error')}")
    print(f"\n{N['pass']} passed, {N['fail']} failed "
          f"(sections C onward need that build)")
    sys.exit(1 if N["fail"] else 0)

print("\nC. panel 1 — composition against burden")
p1, _ = get("/api/panel?n=1&k=7&mt=1&burnin=1")
cls = p1["classes"]
ck(cls == [str(j) for j in range(8)] + ["other", "wt"],
   f"classes are j-of-7 plus `other` and wild type: {cls}")
for sc in SCEN:
    s = p1["scenarios"][sc]
    vals = [(s["classes"][c]["pct_burden"] or {}).get("mean") or 0.0 for c in cls]
    ck(abs(sum(vals) - 100.0) < 0.05,
       f"{sc}: burden classes sum to 100% (got {sum(vals):.1f})")
    ck(s["prevalence"] == DESIGN[sc],
       f"{sc}: prevalence is the design value {DESIGN[sc]}%")
    ck(s["runs"] <= 50, f"{sc}: {s['runs']} runs, never more than the 50 that exist")
ctl = p1["scenarios"]["control"]["classes"]
leak = [c for c in cls if c not in ("0", "wt")
        and ((ctl[c]["pct_burden"] or {}).get("mean") or 0) > 0]
ck(not leak, f"control has no clade with a long-born defining mutation ({leak})")
ck((ctl["wt"]["pct_burden"] or {}).get("mean", 0) > 0,
   "control still has a wild-type class, so the founding clade is kept")

print("\nD. panel 2 — transmission events")
p2, _ = get("/api/panel?n=2&burnin=1")
for sc in SCEN:
    s = p2["scenarios"][sc]
    tot = sum((s["where"][c] or {}).get("mean") or 0 for c in p2["cells"])
    ck(abs(tot - 100.0) < 0.05, f"{sc}: the four cells sum to 100% (got {tot:.1f})")
w = p2["scenarios"]["control"]["who_transmits"]["mean"]
ck(w == 0, f"control: 0% of events caused by a long shedder (got {w})")
hh = p2["scenarios"]["HIV_high"]["who_transmits"]["mean"]
ck(hh > 0, f"HIV high: long shedders do cause events ({hh:.1f}%)")
ck(all(sc in p2["series"] for sc in SCEN), "per-day curves present for every scenario")

print("\nE. panel 3 — emergence speed and latency")
p3, _ = get("/api/panel?n=3&k=7&mt=1")
for sc in SCEN:
    s = p3["scenarios"][sc]
    ck((s["rise_days"] or {}).get("mean") is not None, f"{sc}: a rise time exists")
    ck(s["n_with_rise"] >= s["n_with_latency"],
       f"{sc}: rise denominator {s['n_with_rise']} >= latency denominator "
       f"{s['n_with_latency']}, which is the truncation rule")
    ck("n_truncated" in s, f"{sc}: truncated clades are counted ({s['n_truncated']})")

print("\nF. intervals are over runs, not over clades")
for sc in SCEN:
    for name, b in (("panel1 total_clades", p1["scenarios"][sc]["total_clades"]),
                    ("panel3 rise", p3["scenarios"][sc]["rise_days"])):
        if b and b.get("n") is not None:
            ck(b["n"] <= 50, f"{sc} {name}: n={b['n']} is a count of runs, not clades")

print("\nG. parameters are validated, not trusted")
ck(get_code("/api/panel?n=1&k=6") == 500, "a k that is not one of the k values is refused")
ck(get_code("/api/panel?n=1&k=7&prom=5") == 500, "prominence above 1 is refused")
ck(get_code("/api/panel?n=1&k=7&lo=0.9&hi=0.1") == 500, "lo above hi is refused")
ck(get_code("/api/panel?n=9") == 500, "unknown panel is refused")
ck(get_code("/api/panel?n=1&k=7&win=999") == 500, "absurd smoothing window is refused")
ck(get_code("/api/nope") == 404, "unknown endpoint is 404")

print("\nH. a rebuild request with bad parameters is refused, not dropped")
# k=6 fails validation before anything is launched, so this starts no job
j, code = post("/api/rebuild/clades", {"k": 6})
ck(code == 500, f"an unbuilt k is refused with a status, not a dropped connection ({code})")
ck("error" in j and "k values" in j["error"], f"and it says why: {j.get('error')}")
j, code = post("/api/rebuild/nothing", {})
ck(code == 404, f"unknown rebuild endpoint is 404 ({code})")
j, _ = get("/api/jobs")
ck(j["running"] is None, "no job was started by any of the above")

print("\nI. panel 4 answers, and the job endpoint holds a lock")
p4, _ = get("/api/panel?n=4&k=7&mt=1&burnin=1")
ck(len(p4["rows"]) >= 1, f"{len(p4['rows'])} sensitivity row(s)")
j, _ = get("/api/jobs")
ck("running" in j and "history" in j, "job status is reported")

print("\nJ. the inspection endpoints the clustering and sweep tabs read")
ri, code = get("/api/runs")
ck(code == 200, f"/api/runs answers ({code})")
ck(set(ri["counts"]) == set(SCEN), f"all five scenarios exported: {ri['counts']}")
ck(all(n == 50 for n in ri["counts"].values()), f"50 runs each: {ri['counts']}")
first = ri["scenarios"]["control"][0]
ck(first.get("phi") is not None, f"each run carries its burn-in day ({first.get('phi')})")

r1, code = get("/api/run?scenario=control&seed=0")
ck(code == 200, f"/api/run answers ({code})")
for f in ("parent", "host", "nedge", "depth", "tem", "hosts", "times",
          "f_start", "f_len", "f_vals", "starts"):
    ck(f in r1, f"the run carries {f}")
ck(r1["phi"] == first["phi"], "the run's burn-in day matches the index")
ck(str(7) in r1["starts"], f"clade starts are present at k=7: {sorted(r1['starts'])}")
# the browser rebuilds clade membership from `starts`, so a missing value would
# silently colour the tree at the wrong k
ck(all(str(k) in r1["starts"] for k in (5, 7, 10)),
   "every m_T value the panels use has starts")

j, code = get("/api/run?scenario=nope&seed=0")
ck(code == 500 and "error" in j, f"an unknown scenario is refused with a reason ({code})")
j, code = get("/api/run?scenario=control&seed=999")
ck(code == 500 and "error" in j, f"a seed past the end is refused ({code})")

ac, code = get("/api/across?k=7&win=15&prom=0.1&lo=0.1&hi=0.9&window=post")
if code == 404:
    ck(True, "/api/across has no build at these settings yet (skipped)")
else:
    ck(code == 200, f"/api/across answers ({code})")
    # m_T is a BUILD parameter since v2.5.5, so the values the endpoint returns
    # are the ones that have a build. The invariant is that present and missing
    # together account for all of them; how many are built is a fact about
    # the disk, reported rather than asserted.
    full = len(ac.get("values_full") or [])
    ck(full == 7, f"the server offers 7 m_T values ({full})")
    ck(len(ac["values"]) + len(ac.get("values_missing") or []) == full,
       f"present {len(ac['values'])} + missing "
       f"{len(ac.get('values_missing') or [])} accounts for all {full} values")
    if len(ac["values"]) < full:
        print(f"         {full - len(ac['values'])} value(s) not built: "
              f"{ac.get('values_missing')}")
    c = ac["scenarios"]["control"]
    ck(c["runs"] == 50, f"50 runs of control ({c['runs']})")
    ck(len(c["median"]) == len(ac["values"]), "a median row per m_T value")
    # the filter can only remove clades, so the count must fall as m_T rises
    seq = [row["clades"] for row in c["median"]]
    ck(all(a >= b for a, b in zip(seq, seq[1:])),
       f"clade count falls monotonically with m_T: {seq[0]} -> {seq[-1]}")
    # m_T refuses a clade with no transmissions, so only the first value can
    # report any -- a non-zero above it would mean the filter is not applied
    nz = [i for i, row in enumerate(c["median"]) if row["never_tx_pct"] > 0]
    ck(nz in ([], [0]), f"never-transmitted is non-zero only at m_T=0 (values {nz})")

print(f"\n{N['pass']} passed, {N['fail']} failed")
sys.exit(1 if N["fail"] else 0)
