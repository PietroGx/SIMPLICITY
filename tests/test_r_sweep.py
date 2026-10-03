#!/usr/bin/env python3
'''
Does raising R stop the epidemics fizzling out?

Run #3 lost most of its simulations before the horizon: of 100 seeds per
scenario under the distributional consensus, only 33 (control) to 79
(HIV_high) reached day 1095. Almost none died out and none saturated -- they
stopped on extrande's third condition, `t > 60 and infectious == 0`: a few
hosts still flagged infected, none of them in an infectious compartment, so
transmission was over. At R = 1.03, in a population of 1000, that is the
expected behaviour rather than a failure.

This runs the SAME production scenarios off the SAME frozen table as #3, 10
seeds each, under the distributional consensus only, at R = R_long = 1.05 and
1.10 -- equalising the two so infection duration is the only thing separating
the cohorts. It then prints how often each configuration reaches the horizon,
why the rest stop, and how far they get, against #3's R = 1.03 as the baseline.

The calibrated NSR in the table was fitted at R = 1.03 and is NOT right for
these R values. That is deliberate and harmless: the fade-out rate is a
property of transmission, not of the mutation clock. Read completion rates
from this and nothing else -- no clocks, no clades.

Usage on the HPC, from the repo root:

    python tests/test_r_sweep.py --runner slurm \
        --population-size 5000 --r-values 1.05 1.06 1.07 \
        --slurm-mem 16G

    --population-size N population for the swept runs    (default 1000)
    --i0 N              initial infected; by default scales with the
                        population so starting prevalence matches production
    --r-values A B C    R values to test                 (default 1.05 1.1)
    --seeds N           seeds per scenario per R         (default 10)
    --table-exp-num N   frozen table to reuse            (default 3)
    --exp-num N         experiment number to write under (default 904)
    --analyse-only      skip the runs, re-print the report
    --sequential        one experiment at a time

Runs land in Data/rtest_<R>_N<pop>_<scenario>_#<exp-num>, so sweeps at
different populations never collide. The R=1.03 baseline column comes from run
#3 at population 1000: completion rates stay comparable across populations,
absolute counts do not.

Output is a plain-text report; copy the whole block.
'''
import argparse
import csv
import glob
import os
import statistics
import sys
import threading
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'scripts', 'experiments'))

import pandas as pd

from experiment_script_runner import run_experiment_script
from impact_long_shedders_unbound_config import (
    SETUP_DIR_TEMPLATE, TABLE_FILENAME, build_exp_scenario_settings,
    add_slurm_resource_args, set_slurm_resource_env, prod_exp_name,
    USER_FIXED_PARAMS,
)

CONSENSUS = 'distribution'
BASELINE_R = 1.03
BASELINE_POP = 1000          # run #3's population
HORIZON = 1095.0
SCENARIO_ORDER = ['control', 'SOT', 'HIV_low', 'HIV_high', 'edge_case']


def r_tag(r):
    '''1.05 -> R1p05, so the experiment name is filesystem-safe and readable.'''
    return 'R' + format(float(r), 'g').replace('.', 'p')


def exp_name(r, scenario, pop):
    return f'rtest_{r_tag(r)}_N{pop}_{scenario}'


# ---------------------------------------------------------------- running

def rows_from_table(table_exp_num):
    path = os.path.join(SETUP_DIR_TEMPLATE.format(exp_num=table_exp_num),
                        TABLE_FILENAME)
    if not os.path.isfile(path):
        raise SystemExit(
            f'No frozen table at {path}.\n'
            f'Point --table-exp-num at a completed calibration.')
    df = pd.read_csv(path)
    order = {n: i for i, n in enumerate(SCENARIO_ORDER)}
    rows = [r for _, r in df.iterrows()]
    rows.sort(key=lambda r: order.get(r['scenario_name'], 99))
    return path, rows


def with_r(row, r):
    '''One table row at R = R_long = r. The builder zeroes R_long for control
    (no long shedders), so setting it there is harmless.'''
    out = dict(row)
    out['R'] = float(r)
    if float(out['long_shedders_ratio']) > 0.0:
        out['R_long'] = float(r)
    return out


def with_overrides(builder, extra):
    '''Wrap the pipeline's own settings builder and override a few fixed
    parameters, rather than rebuilding the parameter dict here.'''
    def make_settings():
        varying, fixed, n_seeds = builder()
        fixed = dict(fixed)
        fixed.update(extra)
        return varying, fixed, n_seeds
    return make_settings


def run_all(rows, r_values, exp_num, runner, n_seeds, parallel, overrides, pop):
    jobs = [(r, row) for r in r_values for row in rows]

    def one(r, row):
        name = exp_name(r, row['scenario_name'], pop)
        settings = with_overrides(
            build_exp_scenario_settings(with_r(row, r), n_seeds,
                                        consensus=CONSENSUS), overrides)
        run_experiment_script(runner, exp_num, settings, name)

    if runner != 'slurm' or not parallel or len(jobs) == 1:
        for r, row in jobs:
            one(r, row)
        return []

    # same split as the production dispatcher: every job gets a share of the
    # global cap rather than all of it, and the cap is put back afterwards
    key = 'SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM'
    cap = int(os.environ.get(key, 200))
    budget = max(1, cap // len(jobs))
    print(f'[rtest] submitting {len(jobs)} experiments together: '
          f'{budget} concurrent seeds each, {budget * len(jobs)} of a {cap} cap')

    errors, lock = [], threading.Lock()

    def worker(r, row):
        try:
            one(r, row)
        except BaseException as exc:                      # noqa: BLE001
            with lock:
                errors.append((r, row['scenario_name'], exc))

    previous = os.environ.get(key)
    os.environ[key] = str(budget)
    try:
        threads = [threading.Thread(target=worker, args=(r, row),
                                    name=f"{r_tag(r)}-N{pop}-{row['scenario_name']}")
                   for r, row in jobs]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
    finally:
        if previous is None:
            os.environ.pop(key, None)
        else:
            os.environ[key] = previous
    return errors


# ---------------------------------------------------------------- measuring

def seed_dirs(name, exp_num):
    return sorted(glob.glob(
        f'Data/{name}_#{exp_num}/04_Output/*/seed_*/simulation_trajectory.csv'))


def measure(paths, pop_size):
    '''Completion, stop reason, how far the early ones got, and when phi
    saturates. Mirrors extrande.check_stop_conditions, in its order.

    pop_size must be the population the run actually used: phi is
    (diagnosed + recovered) / size, so a wrong size silently moves every
    burn-in number.'''
    m = {'n': 0, 'full': 0, 'died': 0, 'noinfectious': 0, 'saturated': 0,
         'other': 0, 'stop_days': [], 'peak': [], 'phi_days': [], 'phi_never': 0}
    for p in paths:
        with open(p) as fh:
            rows = list(csv.DictReader(fh))
        if not rows:
            continue
        m['n'] += 1
        m['peak'].append(max(float(r['infected']) for r in rows))

        phi_day = None
        for r in rows:
            phi = min((float(r['diagnosed']) + float(r['recovered'])) / pop_size, 1.0)
            if phi >= 1.0:
                phi_day = float(r['time'])
                break
        if phi_day is None:
            m['phi_never'] += 1
        else:
            m['phi_days'].append(phi_day)

        last = rows[-1]
        t = float(last['time'])
        if t >= HORIZON:
            m['full'] += 1
            continue
        m['stop_days'].append(t)
        infected = float(last['infected'])
        susceptibles = float(last['susceptibles'])
        infectious = float(last['infectious'])
        if infected == 0:
            m['died'] += 1
        elif susceptibles == 0:
            m['saturated'] += 1
        elif infectious == 0 and t > 60.0:
            m['noinfectious'] += 1
        else:
            m['other'] += 1
    return m


def task_seconds(name, exp_num):
    '''Per-simulation wall time, from the .started and .completed signal files
    the Slurm runner drops next to each seeded params file. This is the real
    cost of one task, which is what a longer run scales by -- the sweep's own
    wall clock also carries queue wait and the concurrency cap.'''
    out = []
    pattern = (f'Data/{name}_#{exp_num}/03_Seeded_simulation_parameters/'
               f'*/*.started')
    for started in glob.glob(pattern):
        done = started[:-len('.started')] + '.completed'
        if os.path.exists(done):
            try:
                out.append(os.path.getmtime(done) - os.path.getmtime(started))
            except OSError:
                pass
    return [s for s in out if s >= 0]


def by_r_seconds(r_values, exp_num, pop):
    '''Per-task seconds, kept separate per R. Pooling them hides the thing you
    need: a configuration whose runs die early looks cheap, and the slow
    configuration is the one you would actually commit to.'''
    return {r: [t for s in SCENARIO_ORDER
                for t in task_seconds(exp_name(r, s, pop), exp_num)]
            for r in r_values}


def projection(per_task, cap, seeds, cal1_tasks, cal2_tasks):
    '''Wall-clock estimate for a full pipeline, given one task's cost.

    A stage of N tasks under a cap of C runs in ceil(N/C) waves, each wave
    taking about as long as its slowest member -- so the p95 task is the
    honest multiplier, not the median.
    '''
    import math
    if not per_task:
        return None
    s = sorted(per_task)
    p95 = s[min(len(s) - 1, int(len(s) * 0.95))]
    prod_tasks = 5 * seeds * 2                      # 5 scenarios, 2 consensus
    stages = [('cal_1', cal1_tasks), ('cal_2', cal2_tasks),
              ('production (both consensus)', prod_tasks)]
    rows, total = [], 0.0
    for label, n in stages:
        waves = math.ceil(n / cap) if cap else n
        secs = waves * p95
        total += secs
        rows.append((label, n, waves, secs))
    return {'p95': p95, 'median': statistics.median(s), 'rows': rows,
            'total': total, 'n': len(s)}


def hms(seconds):
    if seconds != seconds:
        return '-'
    seconds = int(round(seconds))
    h, rem = divmod(seconds, 3600)
    m, sec = divmod(rem, 60)
    if h:
        return f'{h}h {m:02d}m'
    if m:
        return f'{m}m {sec:02d}s'
    return f'{sec}s'


def med(xs):
    return statistics.median(xs) if xs else float('nan')


def pct(a, b):
    return 100.0 * a / b if b else float('nan')


def fmt(x, nd=0):
    return '-' if x != x else f'{x:.{nd}f}'


def collect(r_values, exp_num, table_exp_num, pop):
    '''{(R, scenario): metrics}, with the R=1.03 baseline from run #3. The
    baseline ran at BASELINE_POP, the sweep at pop -- phi is normalised by
    whichever applies.'''
    out = {}
    base_prefix = prod_exp_name(CONSENSUS)
    for s in SCENARIO_ORDER:
        paths = sorted(glob.glob(
            f'Data/{base_prefix}_{s}_#{table_exp_num}/04_Output/*/seed_*/'
            f'simulation_trajectory.csv'))
        if paths:
            out[(BASELINE_R, s)] = measure(paths, BASELINE_POP)
    for r in r_values:
        for s in SCENARIO_ORDER:
            paths = seed_dirs(exp_name(r, s, pop), exp_num)
            if paths:
                out[(r, s)] = measure(paths, pop)
    return out


# ---------------------------------------------------------------- report

def report(data, r_values, exp_num, table_exp_num, n_seeds, pop, timing=None):
    L = []
    w = L.append
    rs = [BASELINE_R] + list(r_values)
    present = [r for r in rs if any(k[0] == r for k in data)]

    w('=' * 78)
    w('R SWEEP — does raising R stop the epidemics fizzling out?')
    w('=' * 78)
    w(f'frozen table      : {SETUP_DIR_TEMPLATE.format(exp_num=table_exp_num)}'
      f'/{TABLE_FILENAME}  (calibrated at R={BASELINE_R})')
    w(f'consensus         : {CONSENSUS} only')
    w(f'population        : {pop}'
      + ('' if pop == BASELINE_POP
         else f'   (baseline ran at {BASELINE_POP} — completion is still'
              f' comparable, absolute counts are not)'))
    w(f'R = R_long        : {", ".join(str(r) for r in r_values)}'
      f'   (baseline R={BASELINE_R}, R_long=1.1 from run #{table_exp_num})')
    w(f'seeds per scenario: {n_seeds} swept, baseline as run')
    w('')
    w('NOTE: the NSR in the table was fitted at R=1.03 and is wrong for the')
    w('swept R values. Read completion rates only — no clocks, no clades.')
    w('')

    w('-' * 78)
    w('1. REACHED THE HORIZON (day 1095)      <- the number that matters')
    w('-' * 78)
    head = f'{"scenario":<11}' + ''.join(f'{"R=" + str(r):>16}' for r in present)
    w(head)
    for s in SCENARIO_ORDER:
        line = f'{s:<11}'
        for r in present:
            m = data.get((r, s))
            line += (f'{m["full"]:>7}/{m["n"]:<3}{pct(m["full"], m["n"]):>5.0f}%'
                     if m else f'{"-":>16}')
        w(line)
    line = f'{"ALL":<11}'
    for r in present:
        f_ = sum(data[k]['full'] for k in data if k[0] == r)
        n_ = sum(data[k]['n'] for k in data if k[0] == r)
        line += f'{f_:>7}/{n_:<3}{pct(f_, n_):>5.0f}%'
    w(line)
    w('')

    w('-' * 78)
    w('2. WHY THE REST STOPPED')
    w('-' * 78)
    w(f'{"R":>6}  {"scenario":<11}{"early":>7}{"died out":>10}'
      f'{"no infectious":>15}{"saturated":>11}{"other":>7}{"med stop day":>14}')
    for r in present:
        for s in SCENARIO_ORDER:
            m = data.get((r, s))
            if not m:
                continue
            early = m['n'] - m['full']
            w(f'{r:>6}  {s:<11}{early:>7}{m["died"]:>10}{m["noinfectious"]:>15}'
              f'{m["saturated"]:>11}{m["other"]:>7}{fmt(med(m["stop_days"])):>14}')
        w('')

    w('-' * 78)
    w('3. EPIDEMIC SIZE AND BURN-IN')
    w('-' * 78)
    w(f'{"R":>6}  {"scenario":<11}{"med peak infected":>19}'
      f'{"med day phi=1":>15}{"never saturates":>17}')
    for r in present:
        for s in SCENARIO_ORDER:
            m = data.get((r, s))
            if not m:
                continue
            w(f'{r:>6}  {s:<11}{fmt(med(m["peak"])):>19}'
              f'{fmt(med(m["phi_days"])):>15}'
              f'{m["phi_never"]:>10}/{m["n"]:<6}')
        w('')

    if timing:
        w('-' * 78)
        w('4. TIMING')
        w('-' * 78)
        if timing.get('wall') is not None:
            w(f'sweep wall clock  : {hms(timing["wall"])}'
              f'   ({timing["jobs"]} experiments, {timing["tasks"]} tasks,'
              f' cap {timing["cap"]})')
        proj = timing.get('proj')
        if not proj:
            w('no completed .started/.completed signal pairs found — cannot')
            w('measure per-task time (did this run use the slurm runner?).')
        else:
            w(f'per task, pooled  : median {hms(proj["median"])}, '
              f'p95 {hms(proj["p95"])}   (n={proj["n"]})')
            for r, secs in sorted(timing.get('by_r', {}).items()):
                if secs:
                    ss = sorted(secs)
                    w(f'  R={r:<14}: median {hms(statistics.median(ss))}, '
                      f'p95 {hms(ss[min(len(ss)-1, int(len(ss)*0.95))])}'
                      f'   (n={len(ss)})')
            w('')
            w('Project from the SLOWEST R, not the pooled figure: a')
            w('configuration whose runs die early looks cheap and is not the')
            w('one you would commit to.')
            w('')
            w(f'Projected for a full pipeline at {timing["seeds"]} seeds, '
              f'cap {timing["cap"]}:')
            w(f'  {"stage":<30}{"tasks":>8}{"waves":>8}{"wall":>12}')
            for label, n, waves, secs in proj['rows']:
                w(f'  {label:<30}{n:>8}{waves:>8}{hms(secs):>12}')
            w(f'  {"TOTAL":<30}{"":>8}{"":>8}{hms(proj["total"]):>12}')
            w('')
            w('Waves = ceil(tasks / cap), each costing about the p95 task, so')
            w('this is an upper-ish bound that ignores queue wait. cal_1 and')
            w('cal_2 task counts are run #3 sizes; they also run a shorter')
            w('window than production, so their tasks are cheaper than this')
            w('assumes. Treat the total as an order of magnitude.')
        w('')

    w('-' * 78)
    w('5. VERDICT')
    w('-' * 78)
    base = sum(data[k]['full'] for k in data if k[0] == BASELINE_R)
    basen = sum(data[k]['n'] for k in data if k[0] == BASELINE_R)
    if basen:
        w(f'baseline R={BASELINE_R}: {pct(base, basen):.0f}% reached the horizon')
    for r in r_values:
        f_ = sum(data[k]['full'] for k in data if k[0] == r)
        n_ = sum(data[k]['n'] for k in data if k[0] == r)
        if not n_:
            w(f'         R={r}: no output found')
            continue
        w(f'         R={r}: {pct(f_, n_):.0f}% reached the horizon')
    w('')
    w('A configuration worth calibrating is one where most scenarios clear the')
    w('horizon, so the figures are not built on a biased subset. Watch for the')
    w('opposite failure too: "saturated" above zero means R is now high enough')
    w('to exhaust the susceptibles.')
    w('=' * 78)
    return '\n'.join(L)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--table-exp-num', type=int, default=3)
    p.add_argument('--exp-num', type=int, default=904)
    p.add_argument('--population-size', type=int, default=BASELINE_POP,
                   help='population for the swept runs (default %(default)s)')
    p.add_argument('--i0', type=int, default=None,
                   help='initial infected; default scales the production '
                        'value with the population so starting prevalence '
                        'is unchanged')
    p.add_argument('--seeds', type=int, default=10)
    p.add_argument('--r-values', type=float, nargs='+', default=[1.05, 1.1])
    p.add_argument('--runner', default='slurm',
                   choices=['serial', 'multiprocessing', 'slurm'])
    p.add_argument('--analyse-only', action='store_true')
    p.add_argument('--sequential', action='store_true')
    p.add_argument('--project-seeds', type=int, default=100,
                   help='seeds per scenario to project a full run for')
    p.add_argument('--project-cap', type=int, default=200,
                   help='concurrent tasks assumed when projecting')
    p.add_argument('--project-cal1-tasks', type=int, default=900,
                   help="cal_1 task count (run #3's size)")
    p.add_argument('--project-cal2-tasks', type=int, default=300,
                   help="cal_2 task count (run #3's size)")
    add_slurm_resource_args(p)
    args = p.parse_args()

    table_path, rows = rows_from_table(args.table_exp_num)
    wall = None

    # Scale the initial infected with the population unless told otherwise, so
    # starting prevalence -- and with it the chance of early stochastic death
    # -- is the same as production's 50/1000.
    pop = args.population_size
    base_i0 = int(USER_FIXED_PARAMS['infected_individuals_at_start'])
    i0 = args.i0 if args.i0 is not None else max(
        1, round(base_i0 * pop / BASELINE_POP))
    overrides = {'population_size': pop, 'infected_individuals_at_start': i0}

    if not args.analyse_only:
        set_slurm_resource_env(args.slurm_mem, args.slurm_time)
        print(f'[rtest] table      : {table_path}')
        print(f'[rtest] population : {pop}  (initial infected {i0})')
        if pop > BASELINE_POP and args.slurm_mem == os.environ['SIMPLICITY_SLURM_MEM']:
            # cost scales with the number of infected hosts and the lineages
            # they carry, and 4G was sized on production's 1000
            print(f'[rtest][warn] population is {pop / BASELINE_POP:.0f}x '
                  f'production and --slurm-mem is still the default '
                  f'{args.slurm_mem}. A task killed for memory is a wasted '
                  f'run; consider --slurm-mem 16G --slurm-time 2-00:00:00.')
        print(f'[rtest] R values   : {args.r_values}  (R_long set equal to R)')
        print(f'[rtest] scenarios  : {[r["scenario_name"] for r in rows]}')
        print(f'[rtest] seeds each : {args.seeds}')
        print(f'[rtest] writing to : '
              f'Data/rtest_<R>_N{pop}_<scenario>_#{args.exp_num}')
        t0 = time.monotonic()
        errors = run_all(rows, args.r_values, args.exp_num, args.runner,
                         args.seeds, parallel=not args.sequential,
                         overrides=overrides, pop=pop)
        wall = time.monotonic() - t0
        print(f'[rtest] sweep wall clock: {hms(wall)}')
        for r, s, exc in errors:
            print(f'[rtest][FAILED] R={r} {s}: {type(exc).__name__}: {exc}')

    data = collect(args.r_values, args.exp_num, args.table_exp_num, pop)
    if not data:
        raise SystemExit('No output found to measure. Did the runs complete?')

    by_r = by_r_seconds(args.r_values, args.exp_num, pop)
    per_task = [t for v in by_r.values() for t in v]
    timing = {
        'wall': wall,
        'jobs': len(args.r_values) * len(rows),
        'tasks': len(args.r_values) * len(rows) * args.seeds,
        'cap': args.project_cap,
        'seeds': args.project_seeds,
        'by_r': by_r,
        'proj': projection(per_task, args.project_cap, args.project_seeds,
                           args.project_cal1_tasks, args.project_cal2_tasks),
    }

    text = report(data, args.r_values, args.exp_num, args.table_exp_num,
                  args.seeds, pop, timing)
    print('\n' + text)

    out = os.path.join('Data', f'rtest_report_N{pop}_#{args.exp_num}.txt')
    try:
        os.makedirs('Data', exist_ok=True)
        with open(out, 'w') as fh:
            fh.write(text + '\n')
        print(f'\n[rtest] report also written to {out}')
    except OSError as exc:
        print(f'\n[rtest] could not write the report file: {exc}')


if __name__ == '__main__':
    main()
