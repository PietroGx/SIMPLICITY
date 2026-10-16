#!/usr/bin/env python3
'''What one simulation costs: peak RSS and wall time, measured, not guessed.

    python tests/profile_simulation_cost.py --ladder
    python tests/profile_simulation_cost.py --population-size 5000 \
        --infected-at-start 100 --nsr 1e-3 --final-time 1095

WHY

The Slurm request per task is 4G / 2 days, and both numbers are scar tissue:
memory was raised until OOM kills went 284 -> 72 -> 4 -> 0, and the walltime
was raised from 1 day after edge_case cohorts were killed mid-simulation.
Neither has been re-measured since, and nothing in the repo records what a task
ACTUALLY costs -- so there is no evidence to shrink the request with, only the
memory of the failures that grew it.

This measures it. One simulation, built through the real settings/jobs path and
run through simplicity.runners.unit_run.run_seeded_simulation, in a CHILD
process the parent polls.

  peak RSS  /proc/<pid>/status VmHWM, sampled while it runs -- so a run that is
            killed on the cap or the timeout still reports what it reached.
  rlimit    the child gets RLIMIT_AS, so a runaway dies with MemoryError
            instead of taking the machine into swap. A development box has a
            few GB free; a compute node has more, and the point of the ladder
            is to find the SHAPE of the growth, which does not need the top
            rung to fit.
  reached   the window the run ACTUALLY simulated, plus its infection and
            lineage counts. A row is meaningless without them: final_time is
            a cap, not a promise -- a faster epidemic burns through the
            susceptible pool and ends early, so at R=1.3 a 365-day request
            stops at 37 days and would be reported as a 365-day cost.

Measured here (one seed each), RSS tracks the intra-host LINEAGE count and
almost nothing else -- not population_size, which at R=1.03 does not change
how many infections occur, and not the infection count:

    1,359 lineages  301 MB     25,115 lineages  401 MB
   14,895 lineages  336 MB    113,294 lineages  761 MB

about 293 MB of interpreter and imports plus ~4.2 KB per lineage.

Nothing here is on the pipeline's read path.
'''
import argparse
import json
import os
import resource
import shutil
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)

import simplicity.dir_manager as dm
import simplicity.jobs as jobs
import simplicity.settings_manager as sm

NAME = 'cost_probe'
POLL_S = 0.2


# =============================================================================
# CHILD -- one simulation, under a memory cap
# =============================================================================
def run_child(root, cap_bytes):
    if cap_bytes:
        resource.setrlimit(resource.RLIMIT_AS, (cap_bytes, cap_bytes))
    dm.set_data_dir(root)
    sm._data_dir = root
    from simplicity.runners.unit_run import run_seeded_simulation
    group, record = next(iter(jobs.all_repeats(NAME)))   # the one repeat
    run_seeded_simulation(NAME, group, record['index'])


# =============================================================================
# PARENT -- build, spawn, sample
# =============================================================================
def build(root, params):
    dm.set_data_dir(root)
    sm._data_dir = root
    dm.create_directories(NAME)
    sm.write_experiment_settings(
        NAME, sm.generate_experiment_settings({}, params), 1)
    sm.read_settings_and_write_simulation_parameters(NAME)
    jobs.write_repeats(NAME)


def peak_rss_kb(pid):
    try:
        with open(f'/proc/{pid}/status') as fh:
            for line in fh:
                if line.startswith('VmHWM:'):
                    return int(line.split()[1])
    except (FileNotFoundError, ProcessLookupError, ValueError):
        pass
    return 0


def measure(params, cap_gb, timeout_s):
    '''Returns {rss_mb, seconds, outcome}. Never raises on a failed child:
    a run that does not fit is a measurement, not an error.'''
    root = tempfile.mkdtemp(prefix='cost_probe_')
    try:
        build(root, params)
        cap = int(cap_gb * 1024 ** 3) if cap_gb else 0
        child = subprocess.Popen(
            [sys.executable, os.path.abspath(__file__), '--child', root,
             str(cap)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
            cwd=REPO)
        started, peak = time.time(), 0
        while child.poll() is None:
            peak = max(peak, peak_rss_kb(child.pid))
            if time.time() - started > timeout_s:
                child.kill()
                child.wait()
                return {'rss_mb': peak / 1024, 'seconds': time.time() - started,
                        'outcome': f'TIMEOUT >{timeout_s}s'}
            time.sleep(POLL_S)
        peak = max(peak, peak_rss_kb(child.pid))
        seconds = time.time() - started
        _, err = child.communicate()
        if child.returncode != 0:
            tail = [l for l in err.strip().splitlines() if l.strip()][-1:] or ['']
            kind = 'OUT OF MEMORY' if 'MemoryError' in err else \
                   f'exit {child.returncode}'
            return {'rss_mb': peak / 1024, 'seconds': seconds,
                    'outcome': f'{kind}: {tail[0][:60]}'}
        return {'rss_mb': peak / 1024, 'seconds': seconds, 'outcome': 'ok',
                **what_ran(root, params)}
    finally:
        dm.set_data_dir(root)        # what_ran reads through it, above
        shutil.rmtree(root, ignore_errors=True)


# =============================================================================
def what_ran(root, params):
    '''The window the run reached and what it built there. Asked of the OUTPUT,
    never assumed from the parameters: final_time is a cap, and a run that
    ended early costs what it ran, not what it was asked for.'''
    import simplicity.output_manager as om
    try:
        sod = dm.get_simulation_output_dirs(NAME)[0]
        ssod = dm.get_seeded_simulation_output_dirs(sod)[0]
        return {'reached': om.read_final_time(ssod),
                'asked': params['final_time'],
                'infections': len(om.read_individuals_data(ssod)),
                'lineages': len(om.read_phylogenetic_data(ssod))}
    except Exception:                                           # noqa: BLE001
        return {}


def row(label, params, cap_gb, timeout_s):
    got = measure(params, cap_gb, timeout_s)
    note = got['outcome']
    if 'reached' in got:
        short = '' if got['reached'] >= got['asked'] - 1 else \
                f"  ENDED EARLY (asked {got['asked']:.0f}d)"
        note = (f"{got['reached']:6.0f}d {got['infections']:7d} inf "
                f"{got['lineages']:7d} lin   {got['outcome']}{short}")
    print(f"  {label:30s} {got['rss_mb']:7.0f} MB {got['seconds']:8.1f} s  "
          f"{note}", flush=True)
    return {'label': label, **params, **got}


def ladder(args):
    """Population is the axis the request has to be sized against; the window
    is the axis the walltime has to be. Measured separately, because a run at
    the top of both may not fit on a development machine -- and the shape is
    what sizes the request, not one point."""
    base = {'infected_individuals_at_start': args.infected_at_start,
            'nucleotide_substitution_rate': args.nsr,
            'final_time': args.final_time,
            'R': 1.03, 'IH_virus_emergence_rate': 0.1}
    out = []
    print(f"\npopulation axis   (final_time={args.final_time}, "
          f"NSR={args.nsr:g}, start={args.infected_at_start})")
    print(f"  {'config':34s} {'peak RSS':>11s} {'wall':>11s}   outcome")
    for pop in args.populations:
        out.append(row(f'pop {pop}', {**base, 'population_size': pop},
                       args.mem_cap_gb, args.timeout))

    print(f"\nlong-shedder axis (pop {args.populations[-1]}, "
          f"the intra-host lineage cost)")
    print(f"  {'config':34s} {'peak RSS':>11s} {'wall':>11s}   outcome")
    for ratio, tau in ((0.05, 350.23), (0.2, 350.23)):
        out.append(row(
            f'ratio {ratio} tau {tau:g}',
            {**base, 'population_size': args.populations[-1],
             'long_shedders_ratio': ratio, 'tau_3_long': tau,
             'R_long': 1.1, 'sequence_long_shedders': True,
             'nucleotide_substitution_rate_long': args.nsr_long},
            args.mem_cap_gb, args.timeout))
    return out


def main():
    p = argparse.ArgumentParser(
        description='Peak RSS and wall time of one simulation.')
    p.add_argument('--ladder', action='store_true',
                   help='sweep population and long-shedder load')
    p.add_argument('--populations', type=int, nargs='+',
                   default=[500, 1000, 2000, 5000])
    p.add_argument('--population-size', type=int, default=5000)
    p.add_argument('--infected-at-start', type=int, default=100)
    p.add_argument('--nsr', type=float, default=1e-3)
    p.add_argument('--nsr-long', type=float, default=1e-2)
    p.add_argument('--final-time', type=int, default=365)
    p.add_argument('--long-ratio', type=float, default=0.0)
    p.add_argument('--tau-3-long', type=float, default=350.23)
    p.add_argument('--r', type=float, default=1.03)
    p.add_argument('--r-long', type=float, default=1.1)
    p.add_argument('--mem-cap-gb', type=float, default=2.5,
                   help='RLIMIT_AS for the child (0 disables). Keep it below '
                        'free RAM: the point is to kill a runaway, not swap.')
    p.add_argument('--timeout', type=float, default=1800)
    p.add_argument('--json', type=str, default=None)
    args = p.parse_args()

    print(__doc__.strip().splitlines()[0])
    print(f"interpreter {sys.executable}")
    print(f"cap {args.mem_cap_gb} GB, timeout {args.timeout:g}s")
    started = time.time()

    if args.ladder:
        out = ladder(args)
    else:
        params = {'population_size': args.population_size,
                  'infected_individuals_at_start': args.infected_at_start,
                  'nucleotide_substitution_rate': args.nsr,
                  'final_time': args.final_time,
                  'R': args.r, 'IH_virus_emergence_rate': 0.1}
        if args.long_ratio:
            params.update({'long_shedders_ratio': args.long_ratio,
                           'tau_3_long': args.tau_3_long,
                           'R_long': args.r_long,
                           'sequence_long_shedders': True,
                           'nucleotide_substitution_rate_long': args.nsr_long})
        print(f"\n  {'config':34s} {'peak RSS':>11s} {'wall':>11s}   outcome")
        out = [row(f'pop {args.population_size}', params,
                   args.mem_cap_gb, args.timeout)]

    print(f"\ntotal {time.time() - started:.0f}s")
    if args.json:
        with open(args.json, 'w') as fh:
            json.dump(out, fh, indent=2)
        print(f"written: {args.json}")


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--child':
        run_child(sys.argv[2], int(sys.argv[3]))
    else:
        main()
