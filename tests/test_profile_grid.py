#!/usr/bin/env python3
'''Memory and CPU profile across a grid of settings, on the HPC, 10 seeds a cell.

Answers "what costs most, and what gets worse as we scale" with measurements
rather than argument: every task writes its own RSS curve, container-by-container
byte breakdown, extrande event counters and sampled call stacks, and this script
pools them per cell and prints what moved.

    cd ~/SIMPLICITY && git pull
    python tests/test_profile_grid.py --runner slurm --seeds 10 --exp-num 910

    # one knob at a time
    python tests/test_profile_grid.py --populations 1000 5000 --consensus distribution
    python tests/test_profile_grid.py --final-times 365 1095 --populations 2500

    --analyse-only          re-print the report from what is already on disk
    --tracemalloc           add allocation sites by source line (SLOW, see below)

The grid is population x consensus x horizon, each crossed with the scenarios
in --scenarios. Parameters come from the pipeline's own
build_exp_scenario_settings against the frozen table, so a cell is the
production configuration with one thing changed -- not a parameter dict
rebuilt here, which is how a diagnostic drifts from the stage it diagnoses.

Every cell x scenario is one `_scenario_groups` entry in ONE experiment, so the
whole grid is a single Slurm array (chunked only if it would exceed
--max-array-tasks).

Profiling is switched on by this script, for the tasks, by putting
tests/profiling on PYTHONPATH and setting SIMPLICITY_PROFILE_DIR:
simplicity/runners/slurm.py passes its whole environment through sbatch, and
Python imports sitecustomize at startup from PYTHONPATH. Nothing in simplicity/
or scripts/ is touched.

  WITH --runner slurm this is automatic: each task is a fresh interpreter.
  WITH --runner serial/multiprocessing it is NOT -- those run in this process
  (or a fork of it), where sitecustomize already ran before this script could
  set anything. Export the two variables yourself first; the script says so.

--tracemalloc attributes live bytes to the source line that allocated them,
which is the most direct answer to "what should I change". It roughly doubles
memory and slows every allocation, so those tasks are NOT comparable on time
with the others, and the report keeps them separate. Use it on a small grid.
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
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(REPO, 'scripts', 'experiments'))

import pandas as pd

from experiment_script_runner import run_experiment_script
import simplicity.settings_manager as sm
from impact_long_shedders_unbound_config import (
    SETUP_DIR_TEMPLATE, TABLE_FILENAME, CONSENSUS_MODES,
    build_exp_scenario_settings, add_slurm_resource_args,
    set_slurm_resource_env,
)
# the profiler's own loaders: a second copy of this parsing would be a second
# thing to keep in step with the file format
from profile_report import (load, live, peak_mb, baseline_mb, growth_per_day,
                            CONTAINERS)

PROFILING_DIR = os.path.join(HERE, 'profiling')
GRID_EXP_NAME = 'profile_grid'
SCENARIO_ORDER = ['control', 'SOT', 'HIV_low', 'HIV_high', 'edge_case']

# Whether profiling was already live when this script started. With a
# non-slurm runner the simulation shares this interpreter, whose sitecustomize
# ran before main() could set anything -- so only a pre-set environment works.
_PRESET = os.environ.get('SIMPLICITY_PROFILE_DIR')

W = 100


def rule(char='-'):
    print(char * W)


def heading(title):
    rule()
    print(title)
    rule()


def fmt(value, nd=1, dash='-'):
    if value is None:
        return dash
    return f'{value:.{nd}f}'


def hms(seconds):
    if seconds is None:
        return '-'
    seconds = int(round(seconds))
    if seconds < 60:
        return f'{seconds}s'
    if seconds < 3600:
        return f'{seconds // 60}m {seconds % 60:02d}s'
    return f'{seconds // 3600}h {(seconds % 3600) // 60:02d}m'


def med(values):
    return statistics.median(values) if values else None


def p95(values):
    if not values:
        return None
    ordered = sorted(values)
    return ordered[min(len(ordered) - 1, int(math.ceil(0.95 * len(ordered)) - 1))]


# ---------------------------------------------------------------- building

def rows_from_table(table_exp_num):
    path = os.path.join(SETUP_DIR_TEMPLATE.format(exp_num=table_exp_num),
                        TABLE_FILENAME)
    if not os.path.isfile(path):
        raise SystemExit(f'No frozen table at {path}.')
    frame = pd.read_csv(path)
    order = {name: index for index, name in enumerate(SCENARIO_ORDER)}
    rows = sorted((row for _, row in frame.iterrows()),
                  key=lambda row: order.get(row['scenario_name'], 99))
    return path, rows


def with_r(row, r):
    '''One table row at R = R_long = r, so infection duration is the only
    thing separating the cohorts. The builder zeroes R_long for control.'''
    out = dict(row)
    out['R'] = float(r)
    if float(out['long_shedders_ratio']) > 0.0:
        out['R_long'] = float(r)
    return out


def scenario_lookup(rows):
    '''(long_shedders_ratio, tau_3_long) -> scenario. Unique across the five:
    HIV_low and HIV_high share tau but differ in prevalence. Same key as
    test_config_grid uses, so a cell is identified from the parameters file
    that was actually written rather than from a directory name.'''
    return {(round(float(row['long_shedders_ratio']), 6),
             round(float(row['tau_3_long']), 2)): row['scenario_name']
            for row in rows}


def build_cells(populations, consensus_modes, final_times):
    return [(pop, mode, horizon)
            for pop in populations for mode in consensus_modes
            for horizon in final_times]


def cell_label(cell):
    pop, mode, horizon = cell
    return f'N={pop:<5} {mode[:4]:<5} T={horizon:g}'


def cell_short(cell):
    pop, mode, horizon = cell
    return f'N{pop // 1000}k/{mode[:3]}/{horizon:g}'


def build_groups(rows, cells, n_seeds, r, i0_fraction):
    '''One `_scenario_groups` entry per (cell, scenario).

    Each group is all-scalar, so it is exactly one parameter combination;
    generate_experiment_settings applies group scalars AFTER fixed_params, so
    these win. The consensus mode has to go through
    build_exp_scenario_settings rather than being set here, because it
    selects the phenotype model as well as the distance function.
    '''
    base_fixed, groups = None, []
    for cell in cells:
        pop, mode, horizon = cell
        for row in rows:
            _, fixed, _ = build_exp_scenario_settings(
                with_r(row, r), n_seeds, consensus=mode)()
            if base_fixed is None:
                base_fixed = dict(fixed)
            group = dict(fixed)
            group['population_size'] = pop
            group['infected_individuals_at_start'] = max(1, round(pop * i0_fraction))
            group['final_time'] = float(horizon)
            groups.append(group)
    return groups, (base_fixed or {})


def max_array_size():
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
    '''submit_simulations queues the WHOLE array at once, on hold, so a grid
    past the cluster's array limit is refused outright and the release cap
    never gets a chance to apply.'''
    if max_array_tasks <= 0:
        return [groups]
    per_chunk = max(1, max_array_tasks // max(1, n_seeds))
    return [groups[i:i + per_chunk] for i in range(0, len(groups), per_chunk)]


def part_name(index, total):
    return GRID_EXP_NAME if total == 1 else f'{GRID_EXP_NAME}_part{index + 1}'


def enable_profiling(profile_dir, interval, stack_interval, tracemalloc):
    '''Switch the profiler on for the TASKS, via the environment sbatch
    forwards. Returns the directory the tasks will write into.'''
    os.makedirs(profile_dir, exist_ok=True)
    existing = os.environ.get('PYTHONPATH', '')
    if PROFILING_DIR not in existing.split(os.pathsep):
        os.environ['PYTHONPATH'] = (PROFILING_DIR + os.pathsep + existing
                                    if existing else PROFILING_DIR)
    os.environ['SIMPLICITY_PROFILE_DIR'] = profile_dir
    os.environ['SIMPLICITY_PROFILE_INTERVAL_S'] = str(interval)
    os.environ['SIMPLICITY_PROFILE_STACK_INTERVAL_S'] = str(stack_interval)
    os.environ['SIMPLICITY_PROFILE_STACKS'] = '1'
    os.environ['SIMPLICITY_PROFILE_TRACEMALLOC'] = '1' if tracemalloc else '0'
    return profile_dir


def run_grid(rows, cells, exp_num, runner, n_seeds, max_array_tasks, r,
             i0_fraction):
    groups, base_fixed = build_groups(rows, cells, n_seeds, r, i0_fraction)
    chunks = chunk_groups(groups, n_seeds, max_array_tasks)
    for index, chunk in enumerate(chunks):
        name = part_name(index, len(chunks))

        def make_settings(chunk=chunk):
            return ({'_scenario_groups': chunk}, base_fixed, n_seeds)

        print(f'[profile] --- {name}: {len(chunk)} groups, '
              f'{len(chunk) * n_seeds} tasks ---')
        run_experiment_script(runner, exp_num, make_settings, name)


# ---------------------------------------------------------------- collecting

def experiment_names(exp_num):
    names = [f'{GRID_EXP_NAME}_#{exp_num}']
    names += sorted(os.path.basename(path) for path in
                    glob.glob(f'Data/{GRID_EXP_NAME}_part*_#{exp_num}'))
    return [name for name in names if os.path.isdir(os.path.join('Data', name))]


def seeds_on_disk(exp_num, fallback):
    """The seed count the grid was RUN with, read from the experiment rather
    than from --seeds, which has its own default under --analyse-only."""
    for name in ([GRID_EXP_NAME]
                 + [f'{GRID_EXP_NAME}_part{i}' for i in range(1, 100)]):
        try:
            return int(sm.read_n_seeds_file(f'{name}_#{exp_num}')['n_seeds'])
        except (OSError, KeyError, ValueError):
            continue
    return fallback


def path_key(path):
    """The last two components of a seeded-params path, e.g.
    "<param stem>/seed_0009.json".

    The slurm id map records an ABSOLUTE path, written on the cluster
    (/scratch/.../Data/...), so keying on it only ever matches when the report
    runs on the same machine under the same root. Downloading a run and
    reporting on it locally would silently lose every cell. The tail is unique
    within an experiment -- each group is a distinct parameter combination, so
    each has its own stem -- and it travels.
    """
    return os.path.join(os.path.basename(os.path.dirname(path)),
                        os.path.basename(path))


def task_cells(exp_num, rows):
    '''path_key -> (cell, scenario), from the parameters file each simulation
    actually wrote.'''
    lookup = scenario_lookup(rows)
    out = {}
    for name in experiment_names(exp_num):
        for params_file in glob.glob(
                f'Data/{name}/02_Simulation_parameters/*.json'):
            try:
                with open(params_file) as handle:
                    params = json.load(handle)
            except (OSError, ValueError):
                continue
            key = (round(float(params.get('long_shedders_ratio', 0.0)), 6),
                   round(float(params.get('tau_3_long', 0.0)), 2))
            scenario = lookup.get(key)
            if scenario is None:
                continue
            cell = (int(params['population_size']),
                    str(params.get('consensus', '?')),
                    float(params['final_time']))
            stem = os.path.splitext(os.path.basename(params_file))[0]
            seeded_dir = f'Data/{name}/03_Seeded_simulation_parameters/{stem}'
            for seeded in glob.glob(f'{seeded_dir}/*.json'):
                out[path_key(seeded)] = (cell, scenario)
    return out


def task_wall_seconds(exp_num):
    '''Per-simulation seconds from the .started/.completed signal mtimes.
    Covers tasks whose profile is missing, so completion and cost can be
    reported even where the profiler wrote nothing.'''
    out = {}
    for name in experiment_names(exp_num):
        for started in glob.glob(
                f'Data/{name}/03_Seeded_simulation_parameters/*/*.started'):
            done = started[:-len('.started')] + '.completed'
            if not os.path.exists(done):
                continue
            try:
                delta = os.path.getmtime(done) - os.path.getmtime(started)
            except OSError:
                continue
            if delta >= 0:
                out[path_key(started[:-len('.started')])] = delta
    return out


def signal_counts(exp_num):
    started = completed = failed = total = 0
    for name in experiment_names(exp_num):
        base = f'Data/{name}/03_Seeded_simulation_parameters'
        total += len(glob.glob(f'{base}/*/*.json'))
        started += len(glob.glob(f'{base}/*/*.started'))
        completed += len(glob.glob(f'{base}/*/*.completed'))
        failed += len(glob.glob(f'{base}/*/*.failed'))
    return {'total': total, 'started': started, 'completed': completed,
            'failed': failed}


def attach_cells(tasks, exp_num, rows):
    '''Give every profile its cell and scenario, through the slurm id map that
    slurm.job() writes under the same job/task ids the profile is named by.'''
    mapping = task_cells(exp_num, rows)
    walls = task_wall_seconds(exp_num)
    for task in tasks:
        task['cell'] = None
        task['scenario'] = None
        task['wall_signal'] = None
        found = None
        for name in experiment_names(exp_num):
            candidate = os.path.join('Data', name, 'slurm', 'job_id_mapping',
                                     f"{name}_{task['tag']}.csv")
            if os.path.isfile(candidate):
                found = candidate
                break
        if found is None:
            continue
        try:
            with open(found) as handle:
                params_path = handle.read().strip()
        except OSError:
            continue
        key = path_key(params_path)
        entry = mapping.get(key)
        if entry:
            task['cell'], task['scenario'] = entry
        task['wall_signal'] = walls.get(key)
    return tasks


def task_metrics(task):
    '''One row of numbers per profiled process.'''
    rows = live(task)
    slope, span = growth_per_day(task)
    wall = task['wall_signal'] or task['meta'].get('wall_s')
    cpu = task['meta'].get('cpu_s')
    cpu = cpu if isinstance(cpu, (int, float)) else None
    reached = max((s['sim_time'] for s in rows), default=None)
    reactions = max((s['reactions'] for s in rows if s.get('reactions')),
                    default=None)
    thinning = max((s['thinning'] for s in rows if s.get('thinning')),
                   default=None)
    sizes = {}
    if rows:
        for name in CONTAINERS:
            value = rows[-1].get('mb_' + name)
            if value is not None:
                sizes[name] = value
    return {'tag': task['tag'], 'cell': task['cell'],
            'scenario': task['scenario'],
            'traced': bool(task['meta'].get('tracemalloc')),
            'clean': bool(task['meta'].get('exited_cleanly')),
            'peak_mb': peak_mb(task), 'base_mb': baseline_mb(task),
            'mb_per_day': slope, 'days_span': span, 'days_reached': reached,
            'wall_s': wall, 'cpu_s': cpu,
            'cpu_wall': (cpu / wall) if cpu and wall else None,
            's_per_day': (wall / reached) if wall and reached else None,
            'ms_per_reaction': (1000.0 * wall / reactions)
                               if wall and reactions else None,
            'reactions': reactions, 'thinning': thinning,
            'thinning_pct': (100.0 * thinning / reactions)
                            if reactions and thinning else None,
            'sizes': sizes}


# ---------------------------------------------------------------- reporting

def section_header(exp_num, table_path, seeds, cells, counts, profiled, traced):
    rule('=')
    print('MEMORY AND CPU PROFILE GRID')
    rule('=')
    print(f'frozen table : {table_path}')
    print(f'experiment   : {GRID_EXP_NAME}_#{exp_num}')
    print(f'seeds        : {seeds} per scenario per cell')
    print(f'cells        : {len(cells)}   (population x consensus x horizon)')
    print()
    heading('0. COMPLETENESS          read this before anything else')
    print(f"  simulations  : {counts['completed']} completed, "
          f"{counts['failed']} failed, of {counts['total']} submitted")
    print(f'  profiled     : {profiled} processes wrote a profile')
    if traced:
        print(f'  tracemalloc  : {traced} of them, NOT time-comparable')
    if counts['total'] and profiled < counts['completed']:
        print()
        print('  *** fewer profiles than completed simulations. The profiler is')
        print('  *** enabled through PYTHONPATH, which only reaches a FRESH')
        print('  *** interpreter -- so --runner serial/multiprocessing needs the')
        print('  *** environment exported before launching (see --help).')
    if not profiled:
        print()
        print('  *** no profiles at all: nothing below can be computed.')


def section_cost(metrics, cells):
    heading('1. WHAT EACH CELL COSTS          one row per population x consensus x horizon')
    print(f"{'cell':<22}{'n':>4}{'peak MB':>9}{'p95 MB':>8}{'base':>7}"
          f"{'MB/day':>8}{'wall':>9}{'cpu/wall':>9}{'s/day':>8}"
          f"{'ms/react':>9}{'thin%':>7}")
    by_cell = defaultdict(list)
    for row in metrics:
        if row['cell'] and not row['traced']:
            by_cell[row['cell']].append(row)
    for cell in cells:
        rows = by_cell.get(cell)
        if not rows:
            print(f'{cell_label(cell):<22}{"-":>4}   no profiled simulations')
            continue
        peaks = [r['peak_mb'] for r in rows if r['peak_mb']]
        print(f'{cell_label(cell):<22}{len(rows):>4}'
              f'{fmt(med(peaks)):>9}{fmt(p95(peaks)):>8}'
              f"{fmt(med([r['base_mb'] for r in rows if r['base_mb']])):>7}"
              f"{fmt(med([r['mb_per_day'] for r in rows if r['mb_per_day']]), 3):>8}"
              f"{hms(med([r['wall_s'] for r in rows if r['wall_s']])):>9}"
              f"{fmt(med([r['cpu_wall'] for r in rows if r['cpu_wall']]), 2):>9}"
              f"{fmt(med([r['s_per_day'] for r in rows if r['s_per_day']]), 2):>8}"
              f"{fmt(med([r['ms_per_reaction'] for r in rows if r['ms_per_reaction']]), 2):>9}"
              f"{fmt(med([r['thinning_pct'] for r in rows if r['thinning_pct']])):>7}")
    print()
    print('  peak MB is the kernel high-water mark, the number to compare with')
    print('  --mem. base is RSS before the simulation advances a day: interpreter,')
    print('  matrix-exponential table and the pre-allocated individuals dict, a')
    print('  cost every task pays however short it is.')
    print()
    print('  s/day separates "ran long because it had days to simulate" from')
    print('  "ran long because each day got dearer". ms/react does the same')
    print('  against the event count, which is what extrande actually controls;')
    print('  thin% is the share of events that were rejected thinning steps --')
    print('  work done to advance nothing.')
    print()
    print('  cpu/wall above 1 means threads. With --cpus-per-task=1 that is BLAS')
    print('  oversubscribing the allocation; OMP_NUM_THREADS=1 tests it.')


def section_memory(metrics, cells):
    heading('2. WHERE THE MEMORY IS          measured off the containers, median per cell')
    present = [c for c in cells
               if any(r['cell'] == c and r['sizes'] and not r['traced']
                      for r in metrics)]
    if not present:
        print('  no size samples')
        return
    pooled = defaultdict(dict)
    for cell in present:
        rows = [r for r in metrics if r['cell'] == cell and not r['traced']]
        for name in CONTAINERS:
            values = [r['sizes'][name] for r in rows if name in r['sizes']]
            if values:
                pooled[name][cell] = med(values)
    ranked = sorted(pooled, key=lambda n: -max(pooled[n].values()))
    width = 11
    print(f"{'container':<26}" + ''.join(f'{cell_short(c):>{width}}'
                                         for c in present))
    for name in ranked:
        if max(pooled[name].values()) < 0.05:
            continue
        cells_text = ''.join(f"{fmt(pooled[name].get(c)):>{width}}"
                             for c in present)
        print(f'{name:<26}{cells_text}')
    print()
    print('  Megabytes. Shared objects are charged to every container holding')
    print('  them, so a column sums to more than RSS -- read a row against the')
    print('  peak in section 1, and read ACROSS a row to see what scales.')


def section_scaling(metrics, cells):
    heading('3. WHAT SCALES          the largest population against the smallest')
    populations = sorted({c[0] for c in cells})
    if len(populations) < 2:
        print('  only one population in this grid: nothing to compare.')
        print('  Re-run with --populations 1000 5000 to get this section.')
        return
    low, high = populations[0], populations[-1]
    others = sorted({(c[1], c[2]) for c in cells})
    print(f'{"consensus / horizon":<24}{"metric":<22}'
          f'{"N=" + str(low):>12}{"N=" + str(high):>12}{"ratio":>9}')
    fields = [('peak_mb', 'peak MB', 1), ('base_mb', 'baseline MB', 1),
              ('mb_per_day', 'MB per sim day', 3),
              ('wall_s', 'wall seconds', 0),
              ('s_per_day', 'seconds per day', 2),
              ('ms_per_reaction', 'ms per reaction', 2),
              ('reactions', 'reaction count', 0)]
    for mode, horizon in others:
        shown = False
        for key, label, nd in fields:
            values = {}
            for pop in (low, high):
                rows = [r for r in metrics
                        if r['cell'] == (pop, mode, horizon) and not r['traced']
                        and r.get(key)]
                if rows:
                    values[pop] = med([r[key] for r in rows])
            if len(values) < 2 or not values[low]:
                continue
            ratio = values[high] / values[low]
            tag = f'{mode[:4]} T={horizon:g}' if not shown else ''
            shown = True
            print(f'{tag:<24}{label:<22}{fmt(values[low], nd):>12}'
                  f'{fmt(values[high], nd):>12}{ratio:>8.1f}x')
        if shown:
            print()
    print(f'  Population rose {high / low:.1f}x. A metric rising faster than')
    print('  that is superlinear in population -- those are the ones that')
    print('  decide whether a bigger run is affordable. A metric that barely')
    print('  moves is a fixed cost, and no amount of tuning the model will')
    print('  shift it.')


def section_time(tasks, cells, top=16):
    heading('4. WHERE THE TIME GOES          sampled stacks, % of samples per cell')
    per_cell = defaultdict(lambda: defaultdict(int))
    totals = defaultdict(int)
    for task in tasks:
        cell = task.get('cell')
        if cell is None or task['meta'].get('tracemalloc'):
            continue
        for stack, count in task['stacks'].items():
            leaf = stack.split(';')[-1]
            per_cell[cell][leaf] += count
            totals[cell] += count
    present = [c for c in cells if totals.get(c)]
    if not present:
        print('  no stack samples')
        return
    overall = defaultdict(int)
    for cell in present:
        for leaf, count in per_cell[cell].items():
            overall[leaf] += count
    ranked = sorted(overall, key=lambda leaf: -overall[leaf])[:top]
    width = 11
    print(f"{'leaf function (self time)':<40}"
          + ''.join(f'{cell_short(c):>{width}}' for c in present))
    for leaf in ranked:
        shares = ''.join(
            f'{100.0 * per_cell[c][leaf] / totals[c]:>{width - 1}.1f}%'
            for c in present)
        print(f'{leaf[:38]:<40}{shares}')
    print()
    print('  ' + '  '.join(f'{cell_short(c)}={totals[c]}' for c in present)
          + '  samples')
    print()
    print('  Self time: where the interpreter actually was. Read across a row --')
    print('  a function whose share RISES with population is what to attack')
    print('  first. "import" is interpreter startup, not the model.')


def section_inclusive(tasks, top=14):
    heading('5. WHICH STAGE OWNS THE COST          inclusive time, pooled')
    frames = defaultdict(int)
    total = 0
    for task in tasks:
        if task['meta'].get('tracemalloc'):
            continue
        for stack, count in task['stacks'].items():
            for frame in set(stack.split(';')):
                frames[frame] += count
            total += count
    if not total:
        print('  no stack samples')
        return
    print(f"{'function on stack (it, or something it called)':<60}{'%':>8}")
    shown = 0
    for frame, count in sorted(frames.items(), key=lambda kv: -kv[1]):
        if frame == 'import':
            continue
        print(f'{frame[:58]:<60}{100.0 * count / total:>7.1f}%')
        shown += 1
        if shown == top:
            break
    print()
    print(f'  {total} samples. Inclusive shares do not sum to 100: a frame is')
    print('  counted whenever it is anywhere on the stack. Use this to find the')
    print('  stage, then section 4 to find the line.')


def section_allocations(tasks, top=22):
    traced = [t for t in tasks if t['allocs']]
    if not traced:
        return
    heading('6. ALLOCATION SITES          tracemalloc, live bytes by source line')
    pooled = defaultdict(float)
    counts = defaultdict(int)
    for task in traced:
        for megabytes, count, source in task['allocs']:
            pooled[source] = max(pooled[source], megabytes)
            counts[source] = max(counts[source], count)
    print(f"{'source':<58}{'MB':>10}{'objects':>12}")
    for source, megabytes in sorted(pooled.items(), key=lambda kv: -kv[1])[:top]:
        print(f'{source[:56]:<58}{megabytes:>10.2f}{counts[source]:>12}')
    print()
    print(f'  From {len(traced)} task(s) run with --tracemalloc. This is the one')
    print('  view that names the line responsible, so it is what a fix is')
    print('  written against. tracemalloc doubles memory and slows allocation,')
    print('  so these tasks are excluded from every timing above.')


def section_outliers(metrics, limit=12):
    heading('7. WORST INDIVIDUAL TASKS          what a p95 looks like up close')
    rows = [r for r in metrics if r['peak_mb']]
    if not rows:
        print('  nothing to show')
        return
    print(f"{'tag':<16}{'cell':<22}{'scenario':<12}{'peak MB':>9}"
          f"{'wall':>9}{'days':>7}{'MB/day':>8}  exit")
    for row in sorted(rows, key=lambda r: -(r['peak_mb'] or 0))[:limit]:
        label = cell_label(row['cell']) if row['cell'] else '?'
        print(f"{row['tag'][:14]:<16}{label:<22}"
              f"{(row['scenario'] or '?')[:10]:<12}{fmt(row['peak_mb']):>9}"
              f"{hms(row['wall_s']):>9}{fmt(row['days_reached'], 0):>7}"
              f"{fmt(row['mb_per_day'], 3):>8}  "
              f"{'ok' if row['clean'] else 'KILLED/LIVE'}")
    print()
    print('  A KILLED task is the most informative row here: its last sample is')
    print('  the state the kernel or Slurm stopped it in. Cross-check against')
    print('  sacct -j <jobid> --format=JobID,State,MaxRSS,Elapsed to see which')
    print('  limit it hit.')


def write_task_csv(metrics, path):
    '''Every per-task number, so this can be re-analysed without re-running.'''
    fields = ['tag', 'population', 'consensus', 'final_time', 'scenario',
              'traced', 'clean', 'peak_mb', 'base_mb', 'mb_per_day',
              'days_span', 'days_reached', 'wall_s', 'cpu_s', 'cpu_wall',
              's_per_day', 'ms_per_reaction', 'reactions', 'thinning',
              'thinning_pct'] + ['mb_' + name for name in CONTAINERS]
    with open(path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, restval='')
        writer.writeheader()
        for row in metrics:
            out = {key: row.get(key) for key in fields if key in row}
            if row['cell']:
                out['population'], out['consensus'], out['final_time'] = row['cell']
            for name, value in row['sizes'].items():
                out['mb_' + name] = value
            writer.writerow(out)


# ---------------------------------------------------------------- main

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--table-exp-num', type=int, default=3)
    parser.add_argument('--exp-num', type=int, default=910)
    parser.add_argument('--populations', type=int, nargs='+',
                        default=[1000, 2500, 5000])
    parser.add_argument('--consensus', nargs='+', default=['distribution'],
                        choices=list(CONSENSUS_MODES),
                        help='consensus modes to compare (default '
                             '%(default)s; add argmax to price the two '
                             'pipelines against each other)')
    parser.add_argument('--final-times', type=float, nargs='+',
                        default=[1095.0],
                        help='horizons to compare (default %(default)s)')
    parser.add_argument('--scenarios', nargs='+',
                        default=['control', 'HIV_high'],
                        help='scenarios to include; the two extremes by '
                             'default, since control fades out and HIV_high '
                             'runs to the horizon (default %(default)s). '
                             '"all" uses every scenario in the table.')
    parser.add_argument('--r', type=float, default=1.06)
    parser.add_argument('--i0-fraction', type=float, default=0.03)
    parser.add_argument('--seeds', type=int, default=10)
    parser.add_argument('--runner', default='slurm',
                        choices=['serial', 'multiprocessing', 'slurm'])
    parser.add_argument('--analyse-only', action='store_true')
    parser.add_argument('--max-array-tasks', type=int, default=900)
    parser.add_argument('--profile-interval', type=float, default=15.0,
                        help='seconds between memory samples (default '
                             '%(default)s)')
    parser.add_argument('--stack-interval', type=float, default=0.25,
                        help='seconds between stack samples (default '
                             '%(default)s)')
    parser.add_argument('--tracemalloc', action='store_true',
                        help='also capture allocation sites by source line. '
                             'Roughly doubles memory and slows every '
                             'allocation: use a small grid, and read the '
                             'timings from the untraced tasks.')
    parser.add_argument('--profile-dir', default=None,
                        help='where tasks write profiles (default '
                             'Data/profiles_#<exp-num>, which Data/* already '
                             'gitignores)')
    add_slurm_resource_args(parser)
    args = parser.parse_args()

    table_path, all_rows = rows_from_table(args.table_exp_num)
    if args.scenarios == ['all']:
        rows = all_rows
    else:
        wanted = set(args.scenarios)
        rows = [r for r in all_rows if r['scenario_name'] in wanted]
        missing = wanted - {r['scenario_name'] for r in rows}
        if missing:
            raise SystemExit(f'No such scenario in the table: {sorted(missing)}. '
                             f'Have: {[r["scenario_name"] for r in all_rows]}')
    if not rows:
        raise SystemExit('No scenarios selected.')

    cells = build_cells(args.populations, args.consensus, args.final_times)
    profile_dir = args.profile_dir or os.path.join(
        'Data', f'profiles_#{args.exp_num}')

    if not args.analyse_only:
        if args.runner != 'slurm' and not _PRESET:
            raise SystemExit(
                f'--runner {args.runner} runs the simulation in THIS process, '
                'whose sitecustomize already ran, so setting the environment '
                'now would profile nothing. Export it first:\n\n'
                f'  export PYTHONPATH={PROFILING_DIR}:$PYTHONPATH\n'
                f'  export SIMPLICITY_PROFILE_DIR={profile_dir}\n\n'
                'With --runner slurm each task is a fresh interpreter and this '
                'script sets it for you.')
        enable_profiling(profile_dir, args.profile_interval,
                         args.stack_interval, args.tracemalloc)
        set_slurm_resource_env(args.slurm_mem, args.slurm_time)
        groups = len(cells) * len(rows)
        tasks = groups * args.seeds
        per_chunk = (groups if args.max_array_tasks <= 0
                     else max(1, args.max_array_tasks // max(1, args.seeds)))
        arrays = math.ceil(groups / per_chunk)
        print(f'[profile] table      : {table_path}')
        print(f'[profile] cells      : {len(cells)}  '
              f'({len(args.populations)} populations x '
              f'{len(args.consensus)} consensus x '
              f'{len(args.final_times)} horizons)')
        print(f'[profile] scenarios  : {[r["scenario_name"] for r in rows]}')
        print(f'[profile] R          : {args.r}   '
              f'I0 fraction: {args.i0_fraction}')
        print(f'[profile] groups     : {groups}  (one per cell x scenario, '
              f'all in ONE experiment)')
        print(f'[profile] tasks      : {tasks} total, {arrays} array(s) of at '
              f'most {min(per_chunk, groups) * args.seeds}')
        print(f'[profile] profiles   : {profile_dir}')
        print(f'[profile] sampling   : memory every {args.profile_interval}s, '
              f'stacks every {args.stack_interval}s'
              + (', tracemalloc ON' if args.tracemalloc else ''))
        if args.runner == 'slurm':
            limit = max_array_size()
            per_array = min(per_chunk, groups) * args.seeds
            if limit is None:
                print('[profile][warn] could not read MaxArraySize from '
                      'scontrol; if sbatch refuses the array, lower --seeds.')
            elif per_array > limit:
                raise SystemExit(
                    f'[profile] an array of {per_array} tasks exceeds this '
                    f"cluster's MaxArraySize of {limit}. Lower "
                    f'--max-array-tasks to {limit} or less.')
            else:
                print(f'[profile] MaxArraySize: {limit}, each array fits')
        start = time.monotonic()
        run_grid(rows, cells, args.exp_num, args.runner, args.seeds,
                 args.max_array_tasks, args.r, args.i0_fraction)
        print(f'[profile] wall clock : {hms(time.monotonic() - start)}')

    if not os.path.isdir(profile_dir):
        raise SystemExit(f'No profiles at {profile_dir}. '
                         f'Pass --profile-dir if they went elsewhere.')
    tasks = attach_cells(load(profile_dir), args.exp_num, all_rows)
    metrics = [task_metrics(task) for task in tasks]
    counts = signal_counts(args.exp_num)
    seeds = seeds_on_disk(args.exp_num, args.seeds)

    out_path = os.path.join('Data', f'profile_grid_report_{args.exp_num}.txt')
    csv_path = os.path.join('Data', f'profile_grid_tasks_{args.exp_num}.csv')

    import io
    buffer = io.StringIO()
    stdout = sys.stdout
    sys.stdout = buffer
    try:
        section_header(args.exp_num, table_path, seeds, cells, counts,
                       len(tasks), sum(1 for t in tasks if t['allocs']))
        if tasks:
            section_cost(metrics, cells)
            section_memory(metrics, cells)
            section_scaling(metrics, cells)
            section_time(tasks, cells)
            section_inclusive(tasks)
            section_allocations(tasks)
            section_outliers(metrics)
        rule('=')
    finally:
        sys.stdout = stdout
    text = buffer.getvalue()
    print(text)

    try:
        os.makedirs('Data', exist_ok=True)
        with open(out_path, 'w') as handle:
            handle.write(text)
        write_task_csv(metrics, csv_path)
        print(f'[profile] report written to {out_path}')
        print(f'[profile] per-task numbers written to {csv_path}')
        print(f'[profile] raw samples kept under {profile_dir}')
    except OSError as exc:
        print(f'[profile] could not write the report: {exc}')


if __name__ == '__main__':
    main()
