"""Build clade builds over a grid of parameter sets, one at a time.

The server's job lock refuses a second pass rather than queueing it, so a grid
cannot be driven from the Export page. This runs them sequentially:

    tmux new -s sweep
    python -I tests/studio/sweep.py --procs 2

A set whose directory already carries the current code hash and the full run
count is skipped, so an interrupted sweep resumes where it stopped.

The grid is 3 k x 7 m_T x 4 windows x 2 prominences = 168 builds, which is
about 38 hours sequentially. That is an HPC job, not a local one:
`tests/studio/slurm_precompute.sh` runs it as a Slurm array, one build per
task. Use `--list` to see the grid and `--only N` to run one entry.

`--missing` and `--status` answer "which builds are finished?" for the shell,
so that question has one implementation (`state()`) rather than one here and a
directory count in `pack_cache.sh`.

`--exp-num`/`--exp-name` say which RUN OUTPUT a build reads, and are forwarded
to every compute.py call. `compute.params_key` is the clade parameters plus the
code hash and holds no experiment identity, so experiments are kept apart by
directory instead: `--print-paths` resolves the store, and `out_conflict`
refuses an `--out` belonging to another experiment.
"""
import argparse, itertools, json, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
import compute

K = [5, 7, 10]
WIN = [1.0, 7.0, 15.0, 21.0]
PROM = [0.10, 0.20]
# m_T is a BUILD axis since v2.5.5: rule B designation uses it, so each value is
# its own partition. The values come from compute, which the panels also read.
MT = compute.MT_VALUES
LO, HI, WINDOW = 0.10, 0.90, "post"
DEFAULT_STORE = os.path.join(ROOT, "fig2_studio")


def exp_tag(exp_num, exp_name):
    return f"exp_{exp_name}_{int(exp_num)}"


def is_default_exp(exp_num, exp_name):
    return int(exp_num) == 1 and exp_name == compute.EXP_NAME


def store_dir(store, exp_num=1, exp_name=compute.EXP_NAME):
    """The store holding ONE experiment's grid.

    `compute.params_key` carries no experiment identity, so two experiments
    building into one cache directory take the same directory names and
    `done()` reports the other one's finished build as this one's. Separation
    is therefore by directory: experiment 1 of the default experiment keeps the
    historical `<store>/cache`, anything else gets its own store at
    `<store>/exp_<name>_<num>` -- a store, not a renamed cache, so
    `server.py --store` and `panels.py` read it unchanged.
    """
    store = os.path.expanduser(store)
    if is_default_exp(exp_num, exp_name):
        return store
    return os.path.join(store, exp_tag(exp_num, exp_name))


def cache_dir(store, exp_num=1, exp_name=compute.EXP_NAME):
    return os.path.join(store_dir(store, exp_num, exp_name), "cache")


def out_conflict(out, exp_num, exp_name):
    """Why `out` cannot hold this experiment's builds, or None.

    The refusal is what makes the separation hold: `--out` is also passed by
    hand, `slurm_precompute.sh` takes it positionally, and the historical
    default is experiment 1's. Sharing one directory would serve one
    experiment's builds as another's, and no index.json records which
    experiment a build read.
    """
    parts = os.path.abspath(os.path.expanduser(out)).split(os.sep)
    mine = exp_tag(exp_num, exp_name)
    others = sorted({p for p in parts if p.startswith("exp_")} - {mine})
    if is_default_exp(exp_num, exp_name):
        if others:
            return (f"--out {out} is under {others[0]}, another experiment's "
                    f"store, but this builds {exp_name} #{exp_num}")
        return None
    if mine not in parts:
        return (f"--out {out} has no {mine} component, so {exp_name} "
                f"#{exp_num} would build into another experiment's cache; "
                f"pass --store STORE instead (builds go to "
                f"STORE/{mine}/cache) or an --out under {mine}/")
    return None


def grid():
    """Every parameter set in the grid, in a fixed order.

    The order is the contract with `slurm_precompute.sh`, which turns a Slurm
    array index into one of these rather than restating the axes in bash.
    """
    return [dict(k=k, mt=m, win=w, prom=pr, lo=LO, hi=HI, window=WINDOW)
            for k, m, w, pr in itertools.product(K, MT, WIN, PROM)]


STATES = ("ok", "absent", "no-index", "bad-index", "stale-code", "short-runs")


def state(outdir, want_runs):
    """Why a build is or is not finished. `done()` is this returning "ok".

    The single definition of a finished build. A build that died halfway leaves
    a directory with no index.json or a truncated one, so counting directories
    counts those as present; `pack_cache.sh` asks for this instead.
    """
    if not os.path.isdir(outdir):
        return "absent"
    f = os.path.join(outdir, "index.json")
    if not os.path.isfile(f):
        return "no-index"
    try:
        with open(f) as fh:
            d = json.load(fh)
        n = sum(d.get("counts", {}).values())
    except Exception:
        return "bad-index"
    if d.get("code") != compute.code_hash():
        return "stale-code"
    return "ok" if n >= want_runs else "short-runs"


def done(outdir, want_runs):
    return state(outdir, want_runs) == "ok"


def states(out, want_runs, g=None):
    return [state(os.path.join(out, compute.params_key(p)), want_runs)
            for p in (g or grid())]


def dir_count(out):
    """The OLD measure: directories carrying the current hash, built or not."""
    suffix = "__" + compute.code_hash()
    try:
        return sum(1 for e in os.scandir(out)
                   if e.is_dir() and e.name.endswith(suffix))
    except OSError:
        return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--store", default=DEFAULT_STORE,
                    help="store root. One experiment's grid lives in "
                         "<store>/cache, or <store>/exp_<name>_<num>/cache "
                         "for anything but #1 of the default experiment.")
    ap.add_argument("--out", default=None,
                    help="the cache directory itself. Default: --store "
                         "resolved for --exp-num/--exp-name. An --out that "
                         "belongs to another experiment is refused.")
    ap.add_argument("--exp-num", type=int, default=1,
                    help="which run output a build reads; forwarded to "
                         "compute.py, which used to take its own default of 1 "
                         "whatever the shell printed.")
    ap.add_argument("--exp-name", default=compute.EXP_NAME,
                    help="experiment name, forwarded to compute.py.")
    ap.add_argument("--print-paths", action="store_true",
                    help="print `exp_name exp_num rel store cache` for the "
                         "shell scripts and exit, so the cache location has "
                         "one implementation. `rel` is the store subdirectory "
                         "('.' for the default experiment).")
    ap.add_argument("--procs", type=int, default=2)
    ap.add_argument("--seeds", type=int, default=50)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--list", action="store_true",
                    help="print the grid, one parameter set per line, and exit. "
                         "slurm_precompute.sh indexes this.")
    ap.add_argument("--only", type=int, default=None,
                    help="run just grid entry N (0-based). One Slurm array "
                         "task is one entry.")
    ap.add_argument("--missing", action="store_true",
                    help="print the grid index of every build that is not "
                         "done(), one per line, and exit 0. This is what "
                         "pack_cache.sh and --resubmit read.")
    ap.add_argument("--status", action="store_true",
                    help="print the tally per state, the old directory count "
                         "beside it, and a line per unfinished build.")
    a = ap.parse_args()

    sd = store_dir(a.store, a.exp_num, a.exp_name)
    if a.print_paths:
        print(a.exp_name, a.exp_num,
              os.path.relpath(sd, os.path.expanduser(a.store)),
              sd, os.path.join(sd, "cache"))
        return 0
    if a.out is None:
        a.out = os.path.join(sd, "cache")
    clash = out_conflict(a.out, a.exp_num, a.exp_name)
    if clash:
        print(f"refusing: {clash}", file=sys.stderr)
        return 2
    exp = f"{a.exp_name} #{a.exp_num}"

    want = a.seeds * len(compute.SCEN_ORDER)
    g = grid()
    if a.missing or a.status:
        st = states(a.out, want, g)
        miss = [i for i, s in enumerate(st) if s != "ok"]
        if a.missing:
            for i in miss:
                print(i)
            return 0
        print(f"grid {len(g)} builds, code {compute.code_hash()}, "
              f"want {want} runs each, out {a.out}")
        print(f"  experiment {exp}")
        print("  " + "  ".join(f"{s} {st.count(s)}" for s in STATES
                               if st.count(s)))
        print(f"  directories at this hash {dir_count(a.out)}  "
              f"(the old count; it cannot see a half-built one)")
        print(f"  finished {len(g)-len(miss)} of {len(g)}, "
              f"missing {len(miss)}")
        for i in miss:
            print(f"    {i:4d} {st[i]:10s} {compute.params_key(g[i])}")
        return 0
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

    print(f"experiment {exp} -> {a.out}")
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
                "--exp-num", str(a.exp_num), "--exp-name", a.exp_name,
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
