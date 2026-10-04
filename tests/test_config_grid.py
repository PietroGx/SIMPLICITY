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

A third constraint appeared once N rose: phi = (diagnosed + recovered)/N has to
reach 1 before the burn-in opens, and at N = 5000 seeded with 50 it was not
reaching 1 until day 450-900 of 1095, leaving almost nothing to measure. More
hosts to expose, from the same seed.

So a configuration has to clear three bars at once, and this sweeps the grid to
find where that is, reusing run #3's frozen table unchanged under the
distributional consensus only.

ONE SLURM ARRAY. The whole grid goes in as a single experiment via
settings_manager's `_scenario_groups`, one group per (cell, scenario), exactly
as cal_1 submits its 30 groups as 900 tasks in one array. The previous version
submitted 90 separate arrays from 90 threads; Slurm refused most of them, two
thirds of the grid never ran, and the report printed a decision table anyway.
Hence section 0.

The calibrated NSR in that table was fitted at R = 1.03, N = 1000 and is wrong
for every cell here. Deliberate and harmless: all three bars are properties of
transmission, not of the mutation clock. Configuration only -- no clocks, no
clades.

Usage on the HPC, from the repo root:

    python tests/test_config_grid.py --runner slurm --seeds 50

    --populations A B C     default 1000 2500 5000
    --i0-fractions A B C    default 0.01 0.03 0.05  (production is 50/1000)
    --r-values A B          default 1.06 1.07
    --seeds N               per scenario per cell   (default 10)
    --table-exp-num N       frozen table to reuse   (default 3)
    --exp-num N             write under             (default 906)
    --max-array-tasks N     largest single array      (default 900, 0 = all)
    --analyse-only          skip the run, re-print the report

18 cells x 5 scenarios x N seeds: 900 tasks at 10 seeds, 4500 at 50. The grid
is split into arrays of at most --max-array-tasks (default 900, the size cal_1
submits successfully) and they go in one after another, because
submit_simulations queues a whole array at once -- a grid above the cluster's
array or submit limit is refused outright and the release cap never applies.
Each array's size is checked against MaxArraySize before anything is sent.
'''
import argparse
import csv
import glob
import json
import math
import os
import statistics
import subprocess
import sys
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
GRID_EXP_NAME = 'config_grid'
BASELINE = (1.03, 1000, 50)          # run #3, carried as a reference row


def cell_label(cell):
    r, pop, i0 = cell
    return f'N={pop:<5} I0={i0:<4} R={r}'


# ---------------------------------------------------------------- building

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


def build_grid(r_values, populations, i0_fractions):
    '''(R, population, initial infected), i0 rounded from the fraction.'''
    return [(r, pop, max(1, round(pop * frac)))
            for pop in populations for frac in i0_fractions for r in r_values]


def build_groups(rows, cells, n_seeds):
    '''One `_scenario_groups` entry per (cell, scenario).

    Each group is all-scalar, so it is exactly one parameter combination;
    generate_experiment_settings applies group scalars AFTER fixed_params, so
    these win. Parameters come from the pipeline's own
    build_exp_scenario_settings rather than being rebuilt here -- only
    population and starting infections are added on top.
    '''
    base_fixed, groups = None, []
    for cell in cells:
        r, pop, i0 = cell
        for row in rows:
            _, fixed, _ = build_exp_scenario_settings(
                with_r(row, r), n_seeds, consensus=CONSENSUS)()
            if base_fixed is None:
                base_fixed = dict(fixed)
            group = dict(fixed)
            group['population_size'] = pop
            group['infected_individuals_at_start'] = i0
            groups.append(group)
    return groups, (base_fixed or {})


def max_array_size():
    '''The cluster's MaxArraySize, or None if it cannot be read. A single
    array larger than this is refused outright by sbatch.'''
    try:
        out = subprocess.run(['scontrol', 'show', 'config'],
                             stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                             text=True, timeout=30).stdout
    except (OSError, subprocess.SubprocessError):
        return None
    for line in out.splitlines():
        if line.strip().startswith('MaxArraySize'):
            try:
                return int(line.split('=')[1].strip())
            except (IndexError, ValueError):
                return None
    return None


def chunk_groups(groups, n_seeds, max_array_tasks):
    '''Split the groups so no single array exceeds max_array_tasks.

    submit_simulations puts the WHOLE array in the queue at once, on hold, so
    a grid larger than the cluster's array or submit limit is refused outright
    -- the release cap never gets a chance to apply. Chunks are submitted one
    after another, so only one array is ever queued.
    '''
    if max_array_tasks <= 0:
        return [groups]
    per_chunk = max(1, max_array_tasks // max(1, n_seeds))
    return [groups[i:i + per_chunk]
            for i in range(0, len(groups), per_chunk)]


def part_name(index, total):
    return GRID_EXP_NAME if total == 1 else f'{GRID_EXP_NAME}_part{index + 1}'


def run_grid(rows, cells, exp_num, runner, n_seeds, max_array_tasks):
    groups, base_fixed = build_groups(rows, cells, n_seeds)
    chunks = chunk_groups(groups, n_seeds, max_array_tasks)
    print(f'[grid] submitting {len(chunks)} array(s), '
          f'{len(chunks[0]) * n_seeds} tasks each at most, one after another')
    for i, chunk in enumerate(chunks):
        name = part_name(i, len(chunks))

        def make_settings(chunk=chunk):
            return ({'_scenario_groups': chunk}, base_fixed, n_seeds)

        print(f'[grid] --- {name}: {len(chunk)} groups, '
              f'{len(chunk) * n_seeds} tasks ---')
        run_experiment_script(runner, exp_num, make_settings, name)


# ---------------------------------------------------------------- measuring

def measure(paths, pop_size):
    '''Completion, stop reason, how far the early ones got, and when phi
    saturates. Mirrors extrande.check_stop_conditions, in its order.

    pop_size must be the population the run actually used: phi is
    (diagnosed + recovered)/size, so a wrong size silently moves every
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


def grid_params_files(exp_num):
    '''Every parameter file across the grid's parts. The grid may have been
    submitted as several arrays, each its own experiment.'''
    return glob.glob(f'Data/{GRID_EXP_NAME}_#{exp_num}/'
                     f'02_Simulation_parameters/*.json') + \
           glob.glob(f'Data/{GRID_EXP_NAME}_part*_#{exp_num}/'
                     f'02_Simulation_parameters/*.json')


def merge(a, b):
    '''Combine two metric dicts for the same cell and scenario. Only arises if
    a grid was resubmitted under different chunking.'''
    out = dict(a)
    for k, v in b.items():
        out[k] = out[k] + v if isinstance(v, (int, list)) else v
    return out


def scenario_lookup(rows):
    '''(long_shedders_ratio, tau_3_long) -> scenario. Unique across the five:
    HIV_low and HIV_high share tau but differ in prevalence.'''
    return {(round(float(r['long_shedders_ratio']), 6),
             round(float(r['tau_3_long']), 2)): r['scenario_name']
            for r in rows}


def collect(rows, cells, exp_num, table_exp_num):
    '''{cell: {scenario: metrics}} for the single grid experiment, plus the
    run #3 reference row. Each parameter combination is identified from its
    own written parameters file, not from the directory name.'''
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

    lookup = scenario_lookup(rows)
    for params_file in sorted(grid_params_files(exp_num)):
        try:
            with open(params_file) as fh:
                params = json.load(fh)
        except (OSError, ValueError):
            continue
        key = (round(float(params.get('long_shedders_ratio', 0.0)), 6),
               round(float(params.get('tau_3_long', 0.0)), 2))
        scenario = lookup.get(key)
        if scenario is None:
            continue
        cell = (float(params['R']), int(params['population_size']),
                int(params['infected_individuals_at_start']))
        root = os.path.dirname(os.path.dirname(params_file))
        name = os.path.splitext(os.path.basename(params_file))[0]
        paths = sorted(glob.glob(
            f'{root}/04_Output/{name}/seed_*/simulation_trajectory.csv'))
        if paths:
            per = out.setdefault(cell, {}).get(scenario)
            got = measure(paths, cell[1])
            out[cell][scenario] = merge(per, got) if per else got
    return out


def task_seconds_by_cell(rows, exp_num):
    '''Per-simulation seconds, grouped by cell, from the .started/.completed
    signal mtimes next to each seeded params file.'''
    lookup = scenario_lookup(rows)
    by_cell = {}
    for params_file in sorted(grid_params_files(exp_num)):
        try:
            with open(params_file) as fh:
                params = json.load(fh)
        except (OSError, ValueError):
            continue
        key = (round(float(params.get('long_shedders_ratio', 0.0)), 6),
               round(float(params.get('tau_3_long', 0.0)), 2))
        if lookup.get(key) is None:
            continue
        cell = (float(params['R']), int(params['population_size']),
                int(params['infected_individuals_at_start']))
        root = os.path.dirname(os.path.dirname(params_file))
        name = os.path.splitext(os.path.basename(params_file))[0]
        for started in glob.glob(
                f'{root}/03_Seeded_simulation_parameters/{name}/*.started'):
            done = started[:-len('.started')] + '.completed'
            if os.path.exists(done):
                try:
                    dt = os.path.getmtime(done) - os.path.getmtime(started)
                except OSError:
                    continue
                if dt >= 0:
                    by_cell.setdefault(cell, []).append(dt)
    return by_cell


def failed_count(exp_num):
    return len(glob.glob(f'Data/{GRID_EXP_NAME}_#{exp_num}/'
                         f'03_Seeded_simulation_parameters/*/*.failed')) + \
           len(glob.glob(f'Data/{GRID_EXP_NAME}_part*_#{exp_num}/'
                         f'03_Seeded_simulation_parameters/*/*.failed'))


# ---------------------------------------------------------------- helpers

def pooled(per, key):
    return sum(m[key] for m in per.values())


def pooled_list(per, key):
    return [v for m in per.values() for v in m[key]]


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
    if not n:
        return float('nan')
    p = k / n
    return 100.0 * math.sqrt(max(p * (1 - p), 1e-9) / n)


# ---------------------------------------------------------------- report

def report(data, cells, exp_num, table_exp_num, n_seeds, timing=None):
    L = []
    w = L.append
    expected_per_cell = len(SCENARIO_ORDER) * n_seeds
    ordered = ([BASELINE] if BASELINE in data else []) + \
              [c for c in cells if c in data]

    w('=' * 96)
    w('CONFIGURATION GRID — population x starting infections x R')
    w('=' * 96)
    w(f'frozen table : {SETUP_DIR_TEMPLATE.format(exp_num=table_exp_num)}'
      f'/{TABLE_FILENAME}   (fitted at R=1.03, N=1000)')
    w(f'consensus    : {CONSENSUS} only      R_long set equal to R throughout')
    w(f'seeds        : {n_seeds} per scenario per cell, 5 scenarios '
      f'= {expected_per_cell} per cell')
    w('')

    # ---- section 0: is this report even complete? --------------------------
    w('-' * 96)
    w('0. COMPLETENESS          read this before anything else')
    w('-' * 96)
    complete = partial = absent = 0
    missing_lines = []
    for cell in cells:
        per = data.get(cell, {})
        n = pooled(per, 'n') if per else 0
        if n >= expected_per_cell:
            complete += 1
        elif n:
            partial += 1
            missing_lines.append(f'  PARTIAL  {cell_label(cell):<26}'
                                 f'{n} of {expected_per_cell} simulations, '
                                 f'{len(per)} of 5 scenarios')
        else:
            absent += 1
            missing_lines.append(f'  ABSENT   {cell_label(cell):<26}'
                                 f'no output at all')
    got = sum(pooled(data[c], 'n') for c in cells if c in data)
    want = len(cells) * expected_per_cell
    w(f'cells complete : {complete} of {len(cells)}')
    w(f'simulations    : {got} of {want} ({pct(got, want):.0f}%)')
    nfailed = failed_count(exp_num)
    if nfailed:
        w(f'tasks marked .failed: {nfailed}')
    for line in missing_lines:
        w(line)
    if partial or absent:
        w('')
        w('  *** THE TABLES BELOW ARE INCOMPLETE. A partial cell pools a')
        w('  *** different mix of scenarios from a complete one, so its')
        w('  *** numbers are not comparable. Do not choose a configuration')
        w('  *** from this report until every cell is complete.')
    else:
        w('every cell complete.')
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
        n, full = pooled(per, 'n'), pooled(per, 'full')
        weakest = min((pct(m['full'], m['n']) for m in per.values()),
                      default=float('nan'))
        flag = '  *' if cell == BASELINE else (
            '  !' if n < expected_per_cell and cell != BASELINE else '')
        w(f'{cell_label(cell) + flag:<26}{full:>8}/{n:<4}{pct(full, n):>4.0f}%'
          f'{fmt(weakest):>13}%{pooled(per, "saturated"):>11}'
          f'{fmt(med(pooled_list(per, "phi_days"))):>11}'
          f'{fmt(med(pooled_list(per, "window"))):>9}'
          f'{pooled(per, "phi_never"):>7}/{n:<5}')
    w('')
    w('* run #3, for reference.   ! incomplete cell, not comparable.')
    w('  weakest scen = the worst single scenario\'s completion: an uneven')
    w('  cell reintroduces the survival bias that makes scenarios')
    w('  incomparable.  window = days between phi=1 and the run ending, i.e.')
    w('  what is left to measure once the burn-in opens.')
    w(f'  Completion is +/-{se_pct(expected_per_cell // 2, expected_per_cell):.0f}'
      f' points at {expected_per_cell} simulations — smaller gaps are noise.')
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
        early = pooled(per, 'n') - pooled(per, 'full')
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
            w(f'grid wall clock : {hms(timing["wall"])}   '
              f'({timing["tasks"]} tasks in one array)')
        w(f'{"configuration":<26}{"median task":>14}{"p95 task":>12}'
          f'{"full pipeline at 100 seeds":>30}')
        cap = timing['cap']
        waves = (math.ceil(900 / cap) + math.ceil(300 / cap)
                 + math.ceil(5 * timing['seeds'] * 2 / cap))
        for cell in ordered:
            secs = timing.get('by_cell', {}).get(cell)
            if not secs:
                continue
            ss = sorted(secs)
            p95 = ss[min(len(ss) - 1, int(len(ss) * 0.95))]
            w(f'{cell_label(cell):<26}{hms(statistics.median(ss)):>14}'
              f'{hms(p95):>12}{hms(waves * p95):>30}')
        w('')
        w('Pipeline estimate = (cal_1 + cal_2 + production) waves at the p95')
        w('task, waves = ceil(tasks/cap). Calibration would also move to the')
        w('chosen population, which is why it is costed here.')
        w('')

    w('-' * 96)
    w('5. WHAT TO PICK')
    w('-' * 96)
    w('Rank the cells against the three bars, in this order:')
    w('  1. saturated == 0        a cell that saturates is above the cliff')
    w('  2. weakest scenario high an uneven cell is a biased comparison')
    w('  3. window large          days left to measure after phi=1')
    w('Completion alone is not the criterion: a cell can reach the horizon')
    w('often and still open its burn-in too late to measure anything.')
    w('=' * 96)
    return '\n'.join(L)


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--table-exp-num', type=int, default=3)
    p.add_argument('--exp-num', type=int, default=906)
    p.add_argument('--populations', type=int, nargs='+',
                   default=[1000, 2500, 5000])
    p.add_argument('--i0-fractions', type=float, nargs='+',
                   default=[0.01, 0.03, 0.05])
    p.add_argument('--r-values', type=float, nargs='+', default=[1.06, 1.07])
    p.add_argument('--seeds', type=int, default=10)
    p.add_argument('--runner', default='slurm',
                   choices=['serial', 'multiprocessing', 'slurm'])
    p.add_argument('--analyse-only', action='store_true')
    p.add_argument('--max-array-tasks', type=int, default=900,
                   help='largest single Slurm array to submit; the grid is '
                        'split into that many tasks at a time and submitted '
                        'one array after another (default %(default)s, the '
                        'size cal_1 is known to submit successfully). 0 '
                        'submits the whole grid as one array.')
    p.add_argument('--project-cap', type=int, default=200)
    p.add_argument('--project-seeds', type=int, default=100)
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
        print(f'[grid] groups     : {len(cells) * len(rows)}  '
              f'(one per cell x scenario, all in ONE experiment)')
        n_groups = len(cells) * len(rows)
        per_chunk_groups = (n_groups if args.max_array_tasks <= 0
                            else max(1, args.max_array_tasks // max(1, args.seeds)))
        n_arrays = math.ceil(n_groups / per_chunk_groups)
        per_array = min(per_chunk_groups, n_groups) * args.seeds
        print(f'[grid] tasks      : {tasks} total, '
              f'{n_arrays} array(s) of at most {per_array}, '
              f'submitted one after another')
        print(f'[grid] precision  : {len(rows) * args.seeds} simulations per '
              f'cell, completion +/-'
              f'{se_pct(len(rows) * args.seeds // 2, len(rows) * args.seeds):.0f}'
              f' points')
        if args.runner == 'slurm':
            limit = max_array_size()
            if limit is None:
                print('[grid][warn] could not read MaxArraySize from scontrol; '
                      'if sbatch refuses the array, lower --seeds or split '
                      'the grid with --populations.')
            elif per_array > limit:
                raise SystemExit(
                    f'[grid] an array of {per_array} tasks exceeds this '
                    f'cluster\'s MaxArraySize of {limit}. Lower '
                    f'--max-array-tasks to {limit} or less.')
            else:
                print(f'[grid] MaxArraySize: {limit}, each array fits')
        t0 = time.monotonic()
        run_grid(rows, cells, args.exp_num, args.runner, args.seeds,
                 args.max_array_tasks)
        wall = time.monotonic() - t0
        print(f'[grid] wall clock : {hms(wall)}')

    data = collect(rows, cells, args.exp_num, args.table_exp_num)
    if not data:
        raise SystemExit('No output found to measure. Did the run complete?')

    timing = {'wall': wall,
              'tasks': len(cells) * len(rows) * args.seeds,
              'cap': args.project_cap, 'seeds': args.project_seeds,
              'by_cell': task_seconds_by_cell(rows, args.exp_num)}

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
