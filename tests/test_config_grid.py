#!/usr/bin/env python3
'''
Which (population, starting infections, R) should the next production run use?

Two failure modes bound the usable configuration, and raising R trades one for
the other:

  stochastic fade-out  the epidemic runs out of infectious hosts and extrande
                       stops on `t > 60 and infectious == 0`. At R = 1.03,
                       N = 1000 this took 52% of run #3's simulations.
  saturation           infection outruns diagnosis, the susceptible pool
                       empties and extrande stops on `susceptibles == 0`.
                       At R = 1.10, N = 1000 this took all of them.

A third constraint only showed up once N rose: phi = (diagnosed + recovered)/N
has to reach 1 before the burn-in window opens, and at N = 5000 seeded with 50
it was not reaching 1 until day 450-900 of 1095, leaving almost nothing to
measure. More hosts to expose, from the same seed.

So the configuration has to clear three bars at once -- most runs reach the
horizon, no saturation, and phi saturates early enough to leave a window --
and this sweeps the grid to find where that is. It reuses run #3's frozen
table unchanged, under the distributional consensus only.

The calibrated NSR in that table was fitted at R = 1.03, N = 1000 and is wrong
for every cell here. Deliberate and harmless: all three bars are properties of
transmission, not of the mutation clock. Read this report for configuration
only -- no clocks, no clades.

Usage on the HPC, from the repo root:

    python tests/test_config_grid.py --runner slurm

    --populations A B C     default 1000 2500 5000
    --i0-fractions A B C    default 0.01 0.03 0.05  (production is 50/1000)
    --r-values A B          default 1.06 1.07
    --seeds N               per scenario per cell   (default 10)
    --table-exp-num N       frozen table to reuse   (default 3)
    --exp-num N             write under             (default 905)
    --analyse-only          skip the runs, re-print the report
    --sequential            one experiment at a time

Grid size is populations x fractions x R x 5 scenarios x seeds: 18 cells x 5
scenarios x N seeds, so 900 tasks at 10 seeds and 4500 at 50. The task count,
wave count and a wall-clock estimate print before anything is submitted.

At 10 seeds a cell pools 50 simulations and its completion rate carries about
+/-7 points of standard error -- enough to rank 48% against 90%, not enough to
separate 88% from 90%. At 50 seeds it is 250 simulations and about +/-3.

Runs land in Data/rtest_<R>_N<pop>_I<i0>_<scenario>_#<exp-num>.
'''
import argparse
import csv
import glob
import math
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
)

CONSENSUS = 'distribution'
HORIZON = 1095.0
SCENARIO_ORDER = ['control', 'SOT', 'HIV_low', 'HIV_high', 'edge_case']
# run #3, for the reference row
BASELINE = (1.03, 1000, 50)


def r_tag(r):
    return 'R' + format(float(r), 'g').replace('.', 'p')


def exp_name(cell, scenario):
    r, pop, i0 = cell
    return f'rtest_{r_tag(r)}_N{pop}_I{i0}_{scenario}'


def cell_label(cell):
    r, pop, i0 = cell
    return f'N={pop:<5} I0={i0:<4} R={r}'


# ---------------------------------------------------------------- running

def rows_from_table(table_exp_num):
    path = os.path.join(SETUP_DIR_TEMPLATE.format(exp_num=table_exp_num),
                        TABLE_FILENAME)
    if not os.path.isfile(path):
        raise SystemExit(f'No frozen table at {path}.')
    df = pd.read_csv(path)
    order = {n: i for i, n in enumerate(SCENARIO_ORDER)}
    rows = sorted((r for _, r in df.iterrows()),
                  key=lambda r: order.get(r['scenario_name'], 99))
    return path, rows


def with_r(row, r):
    '''One table row at R = R_long = r, so infection duration is the only
    thing separating the cohorts. The builder zeroes R_long for control.'''
    out = dict(row)
    out['R'] = float(r)
    if float(out['long_shedders_ratio']) > 0.0:
        out['R_long'] = float(r)
    return out


def with_overrides(builder, extra):
    '''Wrap the pipeline's own settings builder rather than rebuilding the
    parameter dict here.'''
    def make_settings():
        varying, fixed, n_seeds = builder()
        fixed = dict(fixed)
        fixed.update(extra)
        return varying, fixed, n_seeds
    return make_settings


def build_grid(r_values, populations, i0_fractions):
    '''(R, population, initial infected), i0 rounded from the fraction.'''
    cells = []
    for pop in populations:
        for frac in i0_fractions:
            i0 = max(1, round(pop * frac))
            for r in r_values:
                cells.append((r, pop, i0))
    return cells


def run_all(rows, cells, exp_num, runner, n_seeds, parallel):
    jobs = [(c, row) for c in cells for row in rows]

    def one(cell, row):
        r, pop, i0 = cell
        settings = with_overrides(
            build_exp_scenario_settings(with_r(row, r), n_seeds,
                                        consensus=CONSENSUS),
            {'population_size': pop, 'infected_individuals_at_start': i0})
        run_experiment_script(runner, exp_num, settings,
                              exp_name(cell, row['scenario_name']))

    if runner != 'slurm' or not parallel or len(jobs) == 1:
        for cell, row in jobs:
            one(cell, row)
        return []

    key = 'SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM'
    cap = int(os.environ.get(key, 200))
    budget = max(1, cap // len(jobs))
    print(f'[grid] submitting {len(jobs)} experiments together: '
          f'{budget} concurrent seeds each, {budget * len(jobs)} of a {cap} cap')

    errors, lock = [], threading.Lock()

    def worker(cell, row):
        try:
            one(cell, row)
        except BaseException as exc:                      # noqa: BLE001
            with lock:
                errors.append((cell, row['scenario_name'], exc))

    previous = os.environ.get(key)
    os.environ[key] = str(budget)
    try:
        threads = [threading.Thread(target=worker, args=(c, row),
                                    name=f"{r_tag(c[0])}-N{c[1]}-I{c[2]}")
                   for c, row in jobs]
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

def measure(paths, pop_size):
    '''Completion, stop reason, how far the early ones got, and when phi
    saturates. Mirrors extrande.check_stop_conditions, in its order.

    pop_size must be the population the run actually used: phi is
    (diagnosed + recovered) / size, so a wrong size silently moves every
    burn-in number without any sign that it has.'''
    m = {'n': 0, 'full': 0, 'died': 0, 'noinfectious': 0, 'saturated': 0,
         'other': 0, 'stop_days': [], 'peak': [], 'phi_days': [],
         'phi_never': 0, 'window': []}
    for p in paths:
        with open(p) as fh:
            rows = list(csv.DictReader(fh))
        if not rows:
            continue
        m['n'] += 1
        m['peak'].append(max(float(r['infected']) for r in rows))

        phi_day = None
        for r in rows:
            if min((float(r['diagnosed']) + float(r['recovered'])) / pop_size,
                   1.0) >= 1.0:
                phi_day = float(r['time'])
                break
        last = rows[-1]
        t = float(last['time'])

        if phi_day is None:
            m['phi_never'] += 1
        else:
            m['phi_days'].append(phi_day)
            # days of measurable run left once the burn-in has opened
            m['window'].append(max(0.0, t - phi_day))

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
    '''Per-simulation wall time from the .started/.completed signal mtimes.'''
    out = []
    for started in glob.glob(f'Data/{name}_#{exp_num}/'
                             f'03_Seeded_simulation_parameters/*/*.started'):
        done = started[:-len('.started')] + '.completed'
        if os.path.exists(done):
            try:
                out.append(os.path.getmtime(done) - os.path.getmtime(started))
            except OSError:
                pass
    return [s for s in out if s >= 0]


def collect(cells, exp_num, table_exp_num):
    '''{cell: {scenario: metrics}}, plus the run #3 reference row.'''
    out = {}
    base_prefix = prod_exp_name(CONSENSUS)
    base = {}
    for s in SCENARIO_ORDER:
        paths = sorted(glob.glob(
            f'Data/{base_prefix}_{s}_#{table_exp_num}/04_Output/*/seed_*/'
            f'simulation_trajectory.csv'))
        if paths:
            base[s] = measure(paths, BASELINE[1])
    if base:
        out[BASELINE] = base
    for cell in cells:
        per = {}
        for s in SCENARIO_ORDER:
            paths = sorted(glob.glob(
                f'Data/{exp_name(cell, s)}_#{exp_num}/04_Output/*/seed_*/'
                f'simulation_trajectory.csv'))
            if paths:
                per[s] = measure(paths, cell[1])
        if per:
            out[cell] = per
    return out


def pooled(per_scenario, key):
    return sum(m[key] for m in per_scenario.values())


def pooled_list(per_scenario, key):
    return [v for m in per_scenario.values() for v in m[key]]


def med(xs):
    return statistics.median(xs) if xs else float('nan')


def pct(a, b):
    return 100.0 * a / b if b else float('nan')


def fmt(x, nd=0):
    return '-' if x != x else f'{x:.{nd}f}'


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


def se_pct(k, n):
    '''Standard error of a completion rate, in points. At 50 seeds a 10-point
    difference is about one standard error -- worth printing so the table is
    not over-read.'''
    if not n:
        return float('nan')
    p = k / n
    return 100.0 * math.sqrt(max(p * (1 - p), 1e-9) / n)


# ---------------------------------------------------------------- report

def report(data, cells, exp_num, table_exp_num, n_seeds, timing=None):
    L = []
    w = L.append
    ordered = ([BASELINE] if BASELINE in data else []) + \
              [c for c in cells if c in data]

    w('=' * 96)
    w('CONFIGURATION GRID — population x starting infections x R')
    w('=' * 96)
    w(f'frozen table : {SETUP_DIR_TEMPLATE.format(exp_num=table_exp_num)}'
      f'/{TABLE_FILENAME}   (fitted at R=1.03, N=1000)')
    w(f'consensus    : {CONSENSUS} only      R_long set equal to R throughout')
    w(f'seeds        : {n_seeds} per scenario per cell, 5 scenarios '
      f'= {5 * n_seeds} per cell')
    w('')
    w('NOTE: the NSR is wrong for every cell here, on purpose. Fade-out,')
    w('saturation and phi are transmission properties. Configuration only.')
    w('')

    w('-' * 96)
    w('1. THE DECISION TABLE          all three bars have to clear at once')
    w('-' * 96)
    w(f'{"configuration":<26}{"reached horizon":>18}{"weakest scen":>14}'
      f'{"saturated":>11}{"phi=1 day":>11}{"window":>9}{"never phi=1":>13}')
    for cell in ordered:
        per = data[cell]
        n = pooled(per, 'n')
        full = pooled(per, 'full')
        sat = pooled(per, 'saturated')
        never = pooled(per, 'phi_never')
        weakest = min((pct(m['full'], m['n']) for m in per.values()),
                      default=float('nan'))
        tag = cell_label(cell) + ('  *' if cell == BASELINE else '')
        w(f'{tag:<26}{full:>8}/{n:<4}{pct(full, n):>4.0f}%'
          f'{fmt(weakest):>13}%{sat:>11}'
          f'{fmt(med(pooled_list(per, "phi_days"))):>11}'
          f'{fmt(med(pooled_list(per, "window"))):>9}'
          f'{never:>7}/{n:<5}')
    w('')
    w('* run #3, for reference.  weakest scen = the worst single scenario\'s')
    w('  completion: an uneven cell reintroduces the survival bias that makes')
    w('  scenarios incomparable.  window = days between phi=1 and the run')
    w('  ending, i.e. what is left to measure after the burn-in opens.')
    base_n = pooled(data[BASELINE], 'n') if BASELINE in data else 0
    cell_n = 5 * n_seeds
    w(f'  Completion is +/-{se_pct(cell_n // 2, cell_n):.0f} points at '
      f'{cell_n} seeds — differences below that are noise.')
    w('')

    w('-' * 96)
    w('2. COMPLETION BY SCENARIO (%)')
    w('-' * 96)
    w(f'{"configuration":<26}' + ''.join(f'{s:>12}' for s in SCENARIO_ORDER))
    for cell in ordered:
        per = data[cell]
        line = f'{cell_label(cell):<26}'
        for s in SCENARIO_ORDER:
            m = per.get(s)
            line += f'{fmt(pct(m["full"], m["n"])):>11}%' if m else f'{"-":>12}'
        w(line)
    w('')

    w('-' * 96)
    w('3. WHY THE REST STOPPED')
    w('-' * 96)
    w(f'{"configuration":<26}{"early":>7}{"died out":>10}{"no infectious":>15}'
      f'{"saturated":>11}{"other":>7}{"med stop day":>14}{"med peak":>10}')
    for cell in ordered:
        per = data[cell]
        n = pooled(per, 'n')
        early = n - pooled(per, 'full')
        w(f'{cell_label(cell):<26}{early:>7}{pooled(per, "died"):>10}'
          f'{pooled(per, "noinfectious"):>15}{pooled(per, "saturated"):>11}'
          f'{pooled(per, "other"):>7}'
          f'{fmt(med(pooled_list(per, "stop_days"))):>14}'
          f'{fmt(med(pooled_list(per, "peak"))):>10}')
    w('')

    if timing:
        w('-' * 96)
        w('4. COST')
        w('-' * 96)
        if timing.get('wall') is not None:
            w(f'sweep wall clock : {hms(timing["wall"])}   '
              f'({timing["jobs"]} experiments, {timing["tasks"]} tasks)')
        w(f'{"configuration":<26}{"median task":>14}{"p95 task":>12}'
          f'{"full pipeline at 100 seeds":>30}')
        for cell in ordered:
            secs = timing.get('by_cell', {}).get(cell)
            if not secs:
                continue
            ss = sorted(secs)
            p95 = ss[min(len(ss) - 1, int(len(ss) * 0.95))]
            cap = timing['cap']
            waves = (math.ceil(900 / cap) + math.ceil(300 / cap)
                     + math.ceil(5 * timing['seeds'] * 2 / cap))
            w(f'{cell_label(cell):<26}{hms(statistics.median(ss)):>14}'
              f'{hms(p95):>12}{hms(waves * p95):>30}')
        w('')
        w('Pipeline estimate = (cal_1 + cal_2 + production) waves at the p95')
        w('task, where waves = ceil(tasks/cap). Calibration would also have to')
        w('move to the chosen population, which is why it is costed here.')
        w('')

    w('-' * 96)
    w('5. WHAT TO PICK')
    w('-' * 96)
    w('Rank the cells yourself against the three bars, in this order:')
    w('  1. saturated == 0        a cell that saturates is above the cliff')
    w('  2. weakest scenario high an uneven cell is a biased comparison')
    w('  3. window large          days left to measure after phi=1')
    w('Completion rate alone is not the criterion; a cell can reach the')
    w('horizon often and still open its burn-in too late to measure anything.')
    w('=' * 96)
    return '\n'.join(L)


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--table-exp-num', type=int, default=3)
    p.add_argument('--exp-num', type=int, default=905)
    p.add_argument('--populations', type=int, nargs='+',
                   default=[1000, 2500, 5000])
    p.add_argument('--i0-fractions', type=float, nargs='+',
                   default=[0.01, 0.03, 0.05])
    p.add_argument('--r-values', type=float, nargs='+', default=[1.06, 1.07])
    p.add_argument('--seeds', type=int, default=10)
    p.add_argument('--runner', default='slurm',
                   choices=['serial', 'multiprocessing', 'slurm'])
    p.add_argument('--analyse-only', action='store_true')
    p.add_argument('--sequential', action='store_true')
    p.add_argument('--project-cap', type=int, default=200)
    p.add_argument('--project-seeds', type=int, default=100)
    p.add_argument('--assume-p95-minutes', type=float, default=25.0,
                   help='p95 task minutes used for the pre-flight estimate '
                        'only (default %(default)s, measured at N=5000 R=1.07)')
    add_slurm_resource_args(p)
    args = p.parse_args()

    table_path, rows = rows_from_table(args.table_exp_num)
    cells = build_grid(args.r_values, args.populations, args.i0_fractions)
    wall = None

    if not args.analyse_only:
        set_slurm_resource_env(args.slurm_mem, args.slurm_time)
        tasks = len(cells) * len(rows) * args.seeds
        print(f'[grid] table      : {table_path}')
        print(f'[grid] cells      : {len(cells)}  '
              f'({len(args.populations)} populations x '
              f'{len(args.i0_fractions)} fractions x {len(args.r_values)} R)')
        for c in cells:
            print(f'[grid]              {cell_label(c)}')
        cap = int(os.environ.get(
            'SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM', 200))
        waves = math.ceil(tasks / cap) if cap else tasks
        est = waves * args.assume_p95_minutes * 60
        print(f'[grid] tasks      : {tasks}  '
              f'({len(rows)} scenarios x {args.seeds} seeds per cell)')
        print(f'[grid] estimate   : {waves} waves at cap {cap}, '
              f'~{hms(est)} at a {args.assume_p95_minutes:.0f}min p95 task')
        print(f'[grid]              (p95 was 20m at N=5000 R=1.06 and 32m at '
              f'R=1.07; smaller populations are cheaper)')
        print(f'[grid] precision  : {5 * args.seeds} simulations per cell, '
              f'completion +/-{se_pct(5 * args.seeds // 2, 5 * args.seeds):.0f} points')
        print(f'[grid] writing to : '
              f'Data/rtest_<R>_N<pop>_I<i0>_<scenario>_#{args.exp_num}')
        t0 = time.monotonic()
        errors = run_all(rows, cells, args.exp_num, args.runner, args.seeds,
                         parallel=not args.sequential)
        wall = time.monotonic() - t0
        print(f'[grid] wall clock : {hms(wall)}')
        for cell, s, exc in errors:
            print(f'[grid][FAILED] {cell_label(cell)} {s}: '
                  f'{type(exc).__name__}: {exc}')

    data = collect(cells, args.exp_num, args.table_exp_num)
    if not data:
        raise SystemExit('No output found to measure. Did the runs complete?')

    by_cell = {c: [t for s in SCENARIO_ORDER
                   for t in task_seconds(exp_name(c, s), args.exp_num)]
               for c in cells}
    timing = {'wall': wall, 'jobs': len(cells) * len(rows),
              'tasks': len(cells) * len(rows) * args.seeds,
              'cap': args.project_cap, 'seeds': args.project_seeds,
              'by_cell': by_cell}

    text = report(data, cells, args.exp_num, args.table_exp_num, args.seeds,
                  timing)
    print('\n' + text)

    out = os.path.join('Data', f'config_grid_report_#{args.exp_num}.txt')
    try:
        os.makedirs('Data', exist_ok=True)
        with open(out, 'w') as fh:
            fh.write(text + '\n')
        print(f'\n[grid] report also written to {out}')
    except OSError as exc:
        print(f'\n[grid] could not write the report file: {exc}')


if __name__ == '__main__':
    main()
