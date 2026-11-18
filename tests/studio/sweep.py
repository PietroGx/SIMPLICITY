"""Build clade builds over a grid of parameter sets, one at a time.

The server's job lock refuses a second pass rather than queueing it, so a grid
cannot be driven from the Export page. This runs them sequentially:

    tmux new -s sweep
    python -I tests/studio/sweep.py --procs 2

A set whose directory already carries the current code hash and the full run
count is skipped, so an interrupted sweep resumes where it stopped.

The grid is 3 k x 16 m_T x 4 windows x 2 prominences = 384 builds, which is
about 38 hours sequentially. That is an HPC job, not a local one:
`tests/studio/slurm_precompute.sh` runs it as a Slurm array, one build per
task. Use `--list` to see the grid and `--only N` to run one entry.
"""
import argparse, itertools, json, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
import compute

K = [5, 7, 10]
WIN = [1.0, 7.0, 15.0, 21.0]
PROM = [0.10, 0.20]
# m_T is a BUILD axis since v2.5.5: rule B designation uses it, so each rung is
# its own partition. The ladder comes from compute, which the panels also read.
MT = compute.MT_LADDER
LO, HI, WINDOW = 0.10, 0.90, "post"


def grid():
    """Every parameter set in the grid, in a fixed order.

    The order is the contract with `slurm_precompute.sh`, which turns a Slurm
    array index into one of these rather than restating the axes in bash.
    """
    return [dict(k=k, mt=m, win=w, prom=pr, lo=LO, hi=HI, window=WINDOW)
            for k, m, w, pr in itertools.product(K, MT, WIN, PROM)]


def done(outdir, want_runs):
    f = os.path.join(outdir, "index.json")
    if not os.path.isfile(f):
        return False
    try:
        d = json.load(open(f))
    except Exception:
        return False
    return (d.get("code") == compute.code_hash()
            and sum(d.get("counts", {}).values()) >= want_runs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.expanduser("~/fig2_studio/cache"))
    ap.add_argument("--procs", type=int, default=2)
    ap.add_argument("--seeds", type=int, default=50)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--list", action="store_true",
                    help="print the grid, one parameter set per line, and exit. "
                         "slurm_precompute.sh indexes this.")
    ap.add_argument("--only", type=int, default=None,
                    help="run just grid entry N (0-based). One Slurm array "
                         "task is one entry.")
    a = ap.parse_args()

    want = a.seeds * len(compute.SCEN_ORDER)
    g = grid()
    if a.list:
        for p in g:
            print(f"{p['k']} {p['mt']:g} {p['win']:g} {p['prom']:g} "
                  f"{p['lo']:g} {p['hi']:g} {p['window']}")
        return 0
    if a.only is not None:
        if not 0 <= a.only < len(g):
            print(f"--only {a.only} outside 0..{len(g)-1}")
            return 2
        g = [g[a.only]]
    todo = [p for p in g
            if not done(os.path.join(a.out, compute.params_key(p)), want)]

    print(f"{len(g)} parameter set(s), {len(g)-len(todo)} already built, "
          f"{len(todo)} to run at {want} runs each")
    for p in g:
        mark = "  " if p in todo else "ok"
        print(f"  {mark} {compute.params_key(p)}")
    if a.dry_run:
        return 0

    t0, fails = time.time(), []
    for i, p in enumerate(todo, 1):
        key = compute.params_key(p)
        print(f"\n=== [{i}/{len(todo)}] {key}   elapsed {(time.time()-t0)/60:.0f} min",
              flush=True)
        argv = [sys.executable, "-I", os.path.join(HERE, "compute.py"),
                "--out", a.out, "--procs", str(a.procs), "--seeds", str(a.seeds),
                "--window", p["window"], "--k", str(p["k"]),
                "--mt", str(p["mt"]),
                "--win", str(p["win"]), "--prom", str(p["prom"]),
                "--lo", str(p["lo"]), "--hi", str(p["hi"])]
        s = time.time()
        rc = subprocess.call(argv, cwd=ROOT)
        print(f"--- {key} rc={rc} in {(time.time()-s)/60:.1f} min", flush=True)
        if rc != 0:
            fails.append((key, rc))

    print(f"\nsweep done in {(time.time()-t0)/60:.0f} min, "
          f"{len(todo)-len(fails)}/{len(todo)} built")
    for key, rc in fails:
        print(f"  FAILED rc={rc}  {key}")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
