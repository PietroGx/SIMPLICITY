#!/usr/bin/env python3
"""clade_preprocess.episodes_of and web/episodes.js episodesOf must agree,
censored episodes included, and the censored flag must land where it should.

    python -I tests/studio/test_censored_episodes.py            synthetic, no data
    python -I tests/studio/test_censored_episodes.py --real [N files, 100]

--real reads N of the cached daily-trajectory memos as plain files.
"""
import json, os, subprocess, sys, tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
for _p in (HERE, os.path.join(ROOT, "tests"), os.path.join(ROOT, "scripts", "nature_plots")):
    sys.path.insert(0, _p)
import compute

WIN, PROM, MT = 15.0, 0.10, 0.01


def synthetic():
    """(id, ys, expected censored peak index or None)."""
    rng = np.random.default_rng(20261120)
    t = np.arange(400.)
    ramp = lambda n0, n1, lo, hi: np.linspace(lo, hi, n1 - n0)
    sweep = np.r_[np.zeros(200), ramp(200, 250, 0, 1), np.ones(150)]
    # the 212L shape: a plateau at 1.0 carrying float noise after smoothing
    noisy = compute.smooth(sweep + rng.normal(0, 1e-16, 400), t, WIN)
    bump = np.r_[np.zeros(100), np.hanning(100), np.zeros(200)]
    return [
        ("sweep to fixation", sweep, 249),
        ("sweep, smoothed with 1e-16 noise", noisy, "first-at-max"),
        ("still climbing at the last sample", np.r_[np.zeros(300), ramp(300, 400, 0, 1)], None),
        ("ordinary bump", bump, None),
        ("rise to 0.6 then hold at 0.58", np.r_[np.zeros(100), ramp(100, 150, 0, .6), np.full(250, .58)], 149),
        ("bump, then a sweep", np.maximum(0.3 * bump, sweep), 249),
        ("plateau that comes back down", np.r_[np.zeros(100), ramp(100, 150, 0, 1), np.ones(100), ramp(250, 300, 1, 0), np.zeros(100)], None),
        ("high from the first sample", np.r_[np.full(200, .9), np.full(200, .95)], None),
    ]


def real(n_files):
    """Trajectories straight from the cached daily memos, read as files: no
    key check, so no re-clustering. The trajectories do not depend on the
    episode rule, so a memo from older code is still the right input."""
    import glob
    store = os.path.join(ROOT, "Data", "figures_studio", "exp_impact_long_shedders_2")
    files = sorted(glob.glob(os.path.join(store, "cache", "_panels", "daily_*.json")))
    files = files[::max(1, len(files) // n_files)][:n_files]
    out = []
    for f in files:
        ser = json.load(open(f))
        t = np.asarray(ser["t"], float)
        j0 = (int(np.searchsorted(t, float(ser["phi"]), side="left"))
              if ser.get("phi") is not None else 0)
        for c in ser["clades"]:
            ys = compute.smooth(np.asarray(c["y"], float)[j0:], t[j0:], WIN)
            out.append((f"{os.path.basename(f)}:{c['name']}", ys, "skip"))
    return out


cases = (real(int(sys.argv[2]) if len(sys.argv) > 2 else 100)
         if sys.argv[1:2] == ["--real"] else synthetic())
py = []
for _, ys, _ in cases:
    _, eps = compute.episodes_of(ys, PROM)
    py.append([[a, b, bool(c)] for a, b, c in eps])

js = r"""
const fs = require("fs"), vm = require("vm");
const ctx = {}; vm.createContext(ctx);
vm.runInContext(fs.readFileSync(process.argv[2], "utf8"), ctx);
const cases = JSON.parse(fs.readFileSync(process.argv[3], "utf8"));
fs.writeFileSync(process.argv[4], JSON.stringify(cases.map(ys =>
  ctx.episodesOf(Float64Array.from(ys), %r).eps.map(e => [e.a, e.b, !!e.cens]))));
""" % PROM
with tempfile.TemporaryDirectory(prefix="test_censored_") as d:
    fi, fo, fj = (os.path.join(d, n) for n in ("test_in.json", "test_out.json", "test_run.js"))
    json.dump([ys.tolist() for _, ys, _ in cases], open(fi, "w"))
    open(fj, "w").write(js)
    subprocess.run(["node", fj, os.path.join(HERE, "web", "episodes.js"), fi, fo], check=True)
    jsr = json.load(open(fo))

bad = 0
for (name, ys, want), p, j in zip(cases, py, jsr):
    if p != j:
        bad += 1
        print(f"  FAIL python and browser differ, {name}: py {p} js {j}")
    cz = [e[1] for e in p if e[2]]
    if len(cz) > 1:
        bad += 1
        print(f"  FAIL more than one censored episode, {name}: {p}")
    if want != "skip":
        if want == "first-at-max":
            want = int(np.argmax(ys >= ys.max() * (1 - 1e-9)))
        got = cz[0] if cz else None
        ok = got == want
        bad += not ok
        print(f"  {'ok  ' if ok else 'FAIL'} {name}: censored peak {got}, expected {want}")
    elif cz and ("Clade_0.2.0" in name or os.environ.get("TEST_VERBOSE")):
        print(f"  censored {name}: {p}  max {ys.max():.3f}")
ncz = sum(1 for p in py for e in p if e[2])
print(f"  {ncz} censored episodes over {len(cases)} trajectories")
print(f"\n{len(cases)} trajectories, {bad} failures")
sys.exit(1 if bad else 0)
