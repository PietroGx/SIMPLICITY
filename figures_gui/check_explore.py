#!/usr/bin/env python3
"""/api/explore must BE Figure 5's measurement at Figure 5's rules.

    python -I figures_gui/check_explore.py [--stamp 2] [--seeds 3]

At the build's own window, prominence and rise band, with every exploration
rule at its Figure 5 value, explore.explore() must return the episodes, the
emergence pick, the height, n_found, tem and trunc that panels.run_daily does,
clade by clade, for every scenario. Then each rule is moved and must change
something -- a knob that moves nothing is a knob wired to nothing. Last, a
measurement fact: no clade is transmitted before it emerges.

Reads the local store and the run output; writes nothing but the daily memos
panels already writes.
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, os.pardir))
sys.path.insert(0, HERE)
os.chdir(ROOT)

import compute                                                    # noqa: E402
import panels                                                     # noqa: E402
import sweep                                                      # noqa: E402
import explore                                                    # noqa: E402

fails = []


def check(label, ok, detail=""):
    if not ok:
        fails.append(label)
    print(f"  {'ok  ' if ok else 'FAIL'} {label}{'  ' + detail if detail and not ok else ''}")


def base_q(scen, seed, win, prom, k=7, mt=0.01):
    return dict(scenario=scen, seed=str(seed), k=str(k), mt=str(mt), win=str(win),
                prom=str(prom), lo="0.1", hi="0.9", window="post", mf="0", burnin="1")


KEYS = ("eps", "main", "h", "n_found", "tem", "trunc")


def same_as_run_daily(store, q):
    a = panels.run_daily(store, q)
    b = explore.explore(store, q)
    bad = []
    da = {c["name"]: c for c in a["clades"]}
    for c in b["clades"]:
        ys = c.get("ys") or []
        if ys and max(range(len(ys)), key=ys.__getitem__) == len(ys) - 1 and c["eps"]:
            continue        # still rising on the last day: counted here, not in the build
        ref = da.get(c["name"])
        if ref is None:
            bad.append(f"{c['name']}: not in run_daily")
            continue
        for key in KEYS:
            mine = c.get(key)
            if key == "eps":        # `pre` is explore's own: the build has no hidden-time maximum
                mine = [{k: v for k, v in e.items() if k != "pre"} for e in mine]
            if ref.get(key) != mine:
                bad.append(f"{c['name']}.{key}: {ref.get(key)!r} != {mine!r}")
    if len(da) != len(b["clades"]):
        bad.append(f"{len(da)} clades against {len(b['clades'])}")
    return bad, b


DAYS_TOL = 0.02


def panel_vs_build(store, win, prom, seeds):
    """explore_panel at Figure 5's rules against the build's own n=3b, every
    run of every scenario: the clade set, and per clade nlong, ndef, tem and
    every episode. Exact, except rise days within DAYS_TOL: the daily record
    stores its axis rounded to 2 decimals (compute.run_daily) where the build
    measures on the raw sampling times, so a day difference can differ by
    0.01. Episode starts, heights and amplitudes -- what hidden time and launch
    are made of -- are compared exactly."""
    q = dict(base_q("control", 0, win, prom), seeds=str(seeds))
    ref = panels.panel3_breakouts(store, panels.clean_params(q))["scenarios"]
    got = explore.explore_panel(store, q)["scenarios"]
    bad, n_cl, n_ep, worst, n_rising = [], 0, 0, 0.0, 0
    for sc in compute.SCEN_ORDER:
        rr = {r["seed"]: r for r in ref[sc]["runs"]}
        for r in got[sc]["runs"]:
            b = rr.get(r["seed"])
            if b is None:
                bad.append(f"{sc}/{r['seed']}: not in the build")
                continue
            mine = {c[4]: c[:4] for c in r["clades"]}
            theirs = {c[4]: c[:4] for c in b["clades"]}
            if set(mine) != set(theirs):
                bad.append(f"{sc}/{r['seed']}: clades {sorted(set(mine) ^ set(theirs))[:3]}")
            for name in set(mine) & set(theirs):
                n_cl += 1
                n_ep += len(theirs[name][3])
                a, b = mine[name], theirs[name]
                # one more episode, the last: a clade still rising on the last
                # day, which explore counts and the build does not
                if a[:3] == b[:3] and len(a[3]) == len(b[3]) + 1 and all(
                        x[:3] == y[:3] for x, y in zip(a[3], b[3])):
                    n_rising += 1
                    a = a[:3] + [a[3][:-1]]
                same = a[:3] == b[:3] and len(a[3]) == len(b[3]) and all(
                    x[:3] == y[:3] for x, y in zip(a[3], b[3]))
                if same and a[3]:
                    worst = max(worst, max(abs(x[3] - y[3]) for x, y in zip(a[3], b[3])))
                if not same:
                    bad.append(f"{sc}/{r['seed']}/{name}: {b} != {a}")
    if worst > DAYS_TOL:
        bad.append(f"rise days differ by up to {worst:.3f} d")
    return bad, n_cl, n_ep, worst, n_rising


def signature(res):
    return [(c["name"], tuple((e["a"], e["b"], e["ia"], e["ib"]) for e in c["eps"]), c["tem"])
            for c in res["clades"]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stamp", default="2")
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--panel-seeds", type=int, default=50,
                    help="runs per scenario in the Figure 5 comparison")
    a = ap.parse_args()
    panels.EXP_NUM, panels.EXP_NAME = a.stamp, compute.EXP_NAME
    store = sweep.store_dir(sweep.DEFAULT_STORE, a.stamp, compute.EXP_NAME)

    print("\n=== at Figure 5's rules, explore IS run_daily")
    n = 0
    for win, prom in ((15, 0.1), (7, 0.2), (1, 0.1), (21, 0.2)):
        for scen in compute.SCEN_ORDER:
            for seed in range(a.seeds):
                q = base_q(scen, seed, win, prom)
                try:
                    bad, res = same_as_run_daily(store, q)
                except FileNotFoundError as exc:
                    print(f"  skip {scen}/{seed} win={win} prom={prom}: {exc}")
                    continue
                n += sum(len(c["eps"]) for c in res["clades"])
                check(f"{scen:9s} seed {seed} win={win:>2} prom={prom}", not bad,
                      "; ".join(bad[:3]))
    check("and that compared real episodes", n > 100, f"{n} episodes")

    # a rule may leave one run alone -- SOT seed 0 has no censored peak -- so
    # each must move SOMETHING across these runs, not every run
    print("\n=== every rule moves something (seed 0 of each scenario, win 15, prom 0.1)")
    refs = {scen: (base_q(scen, 0, 15, 0.1),) for scen in compute.SCEN_ORDER}
    refs = {scen: (q0, signature(explore.explore(store, q0))) for scen, (q0,) in refs.items()}
    for name, val in (("kernel", "gauss"), ("kernel", "savgol"), ("kernel", "median"),
                      ("kernel", "none"), ("space", "log"), ("win", "60"),
                      ("prom_abs", "0.01"), ("sep", "90"), ("width", "30"),
                      ("cens", "0"), ("start", "frac"), ("start", "base"),
                      ("trough", "0.05"), ("anchor", "transmission"), ("gate", "start"), ("dip", "low"), ("promref", "peak")):
        moved = [scen for scen, (q0, ref) in refs.items()
                 if signature(explore.explore(store, dict(q0, **{name: val}))) != ref]
        check(f"{name}={val} changes the measurement in {len(moved)} of {len(refs)} runs",
              bool(moved))

    print("\n=== the width bars drawn are the peaks Min width tests (seed 0, win 15)")
    for width in ("0", "45"):
        bad, n_drop = [], 0
        for scen in compute.SCEN_ORDER:
            for kern in ("box", "gauss"):
                res = explore.explore(store, dict(base_q(scen, 0, 15, 0.1), width=width, kernel=kern))
                for c in res["clades"]:
                    kept = sorted(x[0] for x in c.get("wid", []) if x[5])
                    n_drop += sum(1 for x in c.get("wid", []) if not x[5])
                    if kept != sorted(e["b"] for e in c["eps"] if not e["cens"]):
                        bad.append(f"{scen}/{kern}/{c['name']}")
        check(f"width={width}: kept bars are the episodes' peaks, {n_drop} dropped",
              not bad and (width == "0" or n_drop > 0), "; ".join(bad[:3]))

    print("\n=== at the build's rules, explore_panel IS the clade build (every run)")
    for win, prom in ((15, 0.1), (7, 0.2)):
        bad, n_cl, n_ep, worst, n_rising = panel_vs_build(store, win, prom, a.panel_seeds)
        check(f"win={win:>2} prom={prom}: {n_cl} clades, {n_ep} episodes match the build"
              f" (rise days within {worst:.3f} d; {n_rising} still rising at the end, counted here only)",
              not bad and n_ep > 100, "; ".join(bad[:3]))

    print("\n=== a clade is never transmitted before it emerges")
    for scen in compute.SCEN_ORDER:
        res = explore.explore(store, base_q(scen, 0, 15, 0.1))
        pairs = [(c["t_emerge"], c["t_tx1"]) for c in res["clades"] if c["t_tx1"] is not None]
        early = [p for p in pairs if p[1] < p[0] - 1e-9]
        check(f"{scen:9s} {len(pairs)} transmitted clades, none before emergence",
              bool(pairs) and not early, f"{early[:3]}")

    print(f"\n{'all passed' if not fails else f'{len(fails)} FAILED: {fails}'}")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
