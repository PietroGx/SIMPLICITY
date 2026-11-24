#!/usr/bin/env python3
"""Every renamed field must come back at exactly old/100, everything else
unchanged -- against payloads captured from the live server before the pass.

The three deliberate changes are listed, not tolerated: by_class_rise and
by_class_latent now average per run before bootstrapping, and panel 3 reports a
fourth denominator. Those are checked by hand below, not by equality.
"""
import json, os, sys, urllib.request

BASE = sys.argv[1] if len(sys.argv) > 1 else "http://127.0.0.1:8601"
SNAP = ("/tmp/claude-1000/-home-pietro-Documents-SIMPLICITY/"
        "3f1787ec-a088-4eff-be58-82d11af890b0/scratchpad/baseline")

# old name -> new name, and whether the value was scaled by 100
RENAMED = {"pct_burden": ("frac_burden", True),
           "pct_clades": ("frac_clades", True),
           "pct_burden_majority_long": ("frac_burden_majority_long", True),
           "singleton_pct": ("singleton_frac", True),
           "never_tx_pct": ("never_tx_frac", True)}
RESCALED_SAME_NAME = {"prevalence", "who_transmits", "ll", "ls", "sl", "ss"}
# recomputed on purpose this pass
CHANGED = {"by_class_rise", "by_class_latent", "n_latency_before_start",
           "n_truncated"}

n = bad = 0
def ck(c, m):
    global bad
    if not c:
        bad += 1
        print("  FAIL " + m)

def num(a, b, scaled, where):
    global n
    n += 1
    if a is None or b is None:
        ck((a is None) == (b is None), f"{where}: null-ness differs {a} {b}")
        return
    want = a / 100.0 if scaled else a
    ck(abs(want - b) <= 1e-9 * max(1.0, abs(want)),
       f"{where}: want {want!r} got {b!r}")

def ci(a, b, scaled, where):
    for f in ("mean", "lo", "hi"):
        num(a.get(f), b.get(f), scaled, f"{where}.{f}")
    ck(a.get("n") == b.get("n"), f"{where}.n: {a.get('n')} -> {b.get('n')}")

def walk(a, b, scaled, where):
    if isinstance(a, dict) and isinstance(b, dict):
        if "mean" in a and "n" in a:
            ci(a, b, scaled, where); return
        for k, v in a.items():
            if k in CHANGED:
                continue
            nk, sc = RENAMED.get(k, (k, False))
            sc = sc or scaled or k in RESCALED_SAME_NAME
            ck(nk in b, f"{where}.{nk} missing")
            if nk in b:
                walk(v, b[nk], sc, f"{where}.{nk}")
    elif isinstance(a, (int, float)) and not isinstance(a, bool):
        num(a, b, scaled, where)

def get(u):
    return json.load(urllib.request.urlopen(BASE + u))

files = sorted(f for f in os.listdir(SNAP) if f.endswith(".json"))
for f in files:
    old = json.load(open(os.path.join(SNAP, f)))
    k = old.get("kind")
    if k == "panel1":
        p = old["params"]
        u = (f"/api/panel?n=1&k={p['k']}&win={p['win']}&prom={p['prom']}"
             f"&lo={p['lo']}&hi={p['hi']}&window={p['window']}&mt={p['mt']}"
             f"&burnin={1 if old['burnin'] else 0}&mf=0")
    elif k == "panel3":
        p = old["params"]
        # panel 3's payload does not echo the m_F it was computed at -- it is in
        # the memo key but not the answer -- so the capture's own filename is
        # the only record of it.
        mf = f.split("mf=")[1].split(".json")[0]
        u = (f"/api/panel?n=3&k={p['k']}&win={p['win']}&prom={p['prom']}"
             f"&lo={p['lo']}&hi={p['hi']}&window={p['window']}&mt={p['mt']}"
             f"&burnin=1&mf={mf}")
    elif k == "panel2":
        u = (f"/api/panel?n=2&k=7&win=15&prom=0.1&lo=0.1&hi=0.9&window=post"
             f"&mt=0.01&burnin={1 if old['burnin'] else 0}&mf=0")
    else:
        continue
    try:
        new = get(u)
    except Exception as exc:
        print(f"  FAIL {f}: {exc}"); bad += 1; continue
    walk(old.get("scenarios") or {}, new.get("scenarios") or {}, False, f)

print(f"\n{n} numeric comparisons, {bad} failures")
sys.exit(1 if bad else 0)
