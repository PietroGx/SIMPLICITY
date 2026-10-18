#!/usr/bin/env python3
'''Calibrate AND measure every configuration cell, instead of assuming the NSR
does not matter.

    cd ~/SIMPLICITY && git pull
    python tests/test_calibrated_grid.py --runner slurm --exp-base 9200
    python tests/test_calibrated_grid.py --analyse-only --exp-base 9200

Every configuration grid so far (#905-#908) carried this line:

    NOTE: the NSR is wrong for every cell here, on purpose. Fade-out,
    saturation and phi are transmission properties. Configuration only.

That assumption is no longer safe. Fitness depends on the NSR through genome
divergence, and fitness demonstrably moves the epidemic: comparing v2.4.56 with
v2.4.57 under a pinned PYTHONHASHSEED, 4 of 12 epidemic files differ, controlled
against the same commit on both sides (42/42 identical). The mechanism is open.
So this runs each cell's own cal_1 -> cal_2 -> production, and the report puts
the calibrated rates next to the configuration metrics: if the rates move with
the configuration, the earlier grids cannot be read as written.

Each cell is a full pipeline, so the three stages are strictly sequential within
a cell and the cells are independent. They are run concurrently as separate
processes, with the Slurm release cap split between them.

WHAT IS OVERRIDDEN, AND WHY IT HAS TO BE

The pipeline hardcodes the very configuration this grid is choosing --
the chosen one (impact_long_shedders_config.py:88), R_long = R
(SCENARIOS), a 365-day cal_2 window. run_..._pipeline runs each stage as its own
subprocess, so these are applied by tests/cellconfig/sitecustomize.py, which
Python imports in every one of them. Nothing in scripts/ is modified.

R_long is set equal to R, matching every configuration grid run so far -- the R
axis means nothing if production's long shedders transmit at a fixed 1.1 while
the grid measured them at R.

The cal_2 window defaults to 730 rather than 365: at the configurations under
test phi=1 lands around day 500, so a 365-day cal_2 fits entirely below the
burn-in while production measures above it.
'''
import argparse
import csv
import json
import os
import subprocess
import sys
import threading
import time
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
SCRIPTS = os.path.join(REPO, 'scripts', 'experiments')
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)
sys.path.insert(0, SCRIPTS)

import impact_long_shedders_config as cfg
from impact_long_shedders_config import (
    SCENARIOS, exp_name, add_slurm_resource_args,
)

# The grid's cells are full pipelines, so they used to take experiment numbers
# out of the runs' namespace and RUNS.md had to reserve a range. They get their
# own name instead, and with it numbering from 1. Every path in the config is
# derived from EXP_NAME at call time, so rebinding it here moves the cells'
# directories wholesale -- and the same name goes to the children through
# SIMPLICITY_CELL_EXP_NAME, since this process is the only one cellconfig's
# sitecustomize does NOT run in.
GRID_EXP_NAME = 'ils_grid'
cfg.EXP_NAME = GRID_EXP_NAME
# the configuration metrics are measured by test_config_grid's own function, so
# the two grids cannot drift apart on what "saturated" or "phi" mean
from test_config_grid import measure, hms, fmt, se_pct

CELLCONFIG = os.path.join(HERE, 'cellconfig')
PIPELINE = os.path.join(SCRIPTS, 'run_impact_long_shedders_pipeline.py')
W = 104


def rule(char='-'):
    print(char * W)


def heading(title):
    rule()
    print(title)
    rule()


def build_cells(populations, i0_fractions, r_values):
    return [(pop, frac, r)
            for pop in populations for frac in i0_fractions for r in r_values]


def assign_numbers(cells, exp_base):
    """cell -> experiment number, by position: the first cell is exp_base, the
    second exp_base+1, and so on. --cells selects a subset WITHOUT renumbering,
    so a trial run and the full grid agree on which number is which cell."""
    return {cell: exp_base + index for index, cell in enumerate(cells)}


def cell_is_complete(cell, exp_num, args):
    """A cell counts as done when it has a calibration table AND the full
    production run. A partial production has to be redone: pooling it with
    complete cells is the survival bias the grid exists to avoid."""
    if not os.path.isfile(table_path(exp_num)):
        return False
    expected = len(SCENARIOS) * args.exp_seeds * len(args.consensus)
    return len(production_paths_any(exp_num, args.consensus)) >= expected


def cell_label(cell):
    pop, frac, r = cell
    return f'N={pop:<5} I0={max(1, round(pop * frac)):<4} ({frac:.0%}) R={r}'


def cell_short(cell):
    pop, frac, r = cell
    return f'N{pop // 1000}k/{frac:.0%}/{r}'


# ------------------------------------------------------------------- running

def cell_env(cell, args, cap_per_cell):
    pop, frac, r = cell
    env = dict(os.environ)
    path = env.get('PYTHONPATH', '')
    env['PYTHONPATH'] = CELLCONFIG + (os.pathsep + path if path else '')
    env['SIMPLICITY_CELL_POPULATION'] = str(pop)
    env['SIMPLICITY_CELL_I0'] = str(max(1, round(pop * frac)))
    env['SIMPLICITY_CELL_RLONG'] = str(r)
    env['SIMPLICITY_CELL_CAL2_TIME'] = str(args.cal2_final_time)
    env['SIMPLICITY_CELL_NSR_STEPS'] = str(args.nsr_steps)
    env['SIMPLICITY_CELL_EXP_NAME'] = GRID_EXP_NAME
    # the cells share one cluster: without this each would release the full cap
    env['SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM'] = str(cap_per_cell)
    # this script captures the child's stdout to a FILE, which makes Python
    # block-buffer it: the per-cell log sat empty for a whole stage while the
    # cell was running normally. The pipeline's own log is flushed per line and
    # stays the authoritative one, but this should not lag either.
    env['PYTHONUNBUFFERED'] = '1'
    if args.rerun:
        # read by slurm.resuming() and output_manager.setup_output_directory,
        # both several subprocesses down from here
        env['SIMPLICITY_RESUME'] = '1'
    return env


def cell_command(cell, exp_num, args):
    _, _, r = cell
    command = [sys.executable, PIPELINE,
               '--exp-num', str(exp_num),
               '--runner', args.runner,
               '--cal-seeds', str(args.cal_seeds),
               '--exp-seeds', str(args.exp_seeds),
               '--consensus', *args.consensus,
               '--r-cal1', str(r),
               '--r-cal2', str(r),
               '--target-osr-std', str(args.target_osr_std),
               '--target-osr-long', str(args.target_osr_long),
               '--no-compress']
    if args.slurm_mem:
        command += ['--slurm-mem', args.slurm_mem]
    if args.slurm_time:
        command += ['--slurm-time', args.slurm_time]
    return command


def run_cells(cells, numbers, args):
    cap_per_cell = max(1, args.cap // max(1, min(args.concurrency, len(cells))))
    log_dir = os.path.join('Data', 'pipeline_logs')
    os.makedirs(log_dir, exist_ok=True)
    gate = threading.Semaphore(args.concurrency)
    outcomes = {}
    lock = threading.Lock()

    def worker(cell):
        exp_num = numbers[cell]
        with gate:
            command = cell_command(cell, exp_num, args)
            log_path = os.path.join(log_dir, f'calibrated_grid_#{exp_num}.log')
            start = time.monotonic()
            with open(log_path, 'w') as handle:
                handle.write(' '.join(command) + '\n\n')
                handle.flush()
                code = subprocess.call(command, cwd=REPO,
                                       env=cell_env(cell, args, cap_per_cell),
                                       stdout=handle,
                                       stderr=subprocess.STDOUT)
            with lock:
                outcomes[cell] = (code, time.monotonic() - start, exp_num)
            print(f'  [{"ok " if code == 0 else "FAIL"}] {cell_label(cell)}  '
                  f'#{exp_num}  {hms(time.monotonic() - start)}  -> {log_path}',
                  flush=True)

    threads = [threading.Thread(target=worker, args=(c,), daemon=False)
               for c in cells]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    return outcomes


# ----------------------------------------------------------------- collecting

def table_path(exp_num):
    return cfg.table_path(exp_num)


def read_table(exp_num):
    path = table_path(exp_num)
    if not os.path.isfile(path):
        return None
    with open(path) as handle:
        return list(csv.DictReader(handle))


def production_paths(exp_num, consensus):
    root = os.path.join(cfg.run_dir(exp_num, consensus), '04_Output')
    found = []
    for dirpath, _, files in os.walk(root):
        if 'simulation_trajectory.csv' in files:
            found.append(os.path.join(dirpath, 'simulation_trajectory.csv'))
    return found


def production_paths_any(exp_num, consensus_modes):
    """Production lives under one experiment per consensus mode; the scenarios
    are split across sibling directories named by scenario."""
    paths = []
    for consensus in consensus_modes:
        paths += production_paths(exp_num, consensus)
        for scenario in SCENARIOS:
            root = os.path.join(
                'Data',
                f'{exp_name(consensus)}_{scenario["name"]}_#{exp_num}',
                '04_Output')
            for dirpath, _, files in os.walk(root):
                if 'simulation_trajectory.csv' in files:
                    paths.append(os.path.join(dirpath,
                                              'simulation_trajectory.csv'))
    return sorted(set(paths))


def collect(cells, numbers, args):
    data = {}
    for cell in cells:
        exp_num = numbers[cell]
        pop = cell[0]
        rows = read_table(exp_num)
        paths = production_paths_any(exp_num, args.consensus)
        data[cell] = {
            'exp_num': exp_num,
            'table': rows,
            'n_production': len(paths),
            'metrics': measure(paths, pop) if paths else None,
        }
    return data


# ------------------------------------------------------------------ reporting

def section_completeness(cells, data, outcomes, args):
    heading('0. COMPLETENESS          read this before anything else')
    done = sum(1 for c in cells if data[c]['table'])
    print(f'  cells                 : {len(cells)}')
    print(f'  calibration tables    : {done} of {len(cells)}')
    runs = sum(1 for c in cells if data[c]['n_production'])
    print(f'  cells with production : {runs} of {len(cells)}')
    if outcomes:
        failed = [c for c in cells if outcomes.get(c, (1,))[0] != 0]
        print(f'  pipelines exiting non-zero : {len(failed)}')
        for cell in failed:
            print(f'      {cell_label(cell)}  #{data[cell]["exp_num"]}  '
                  f'see Data/pipeline_logs/')
    print()
    print(f'  seeds: {args.cal_seeds} per calibration grid point, '
          f'{args.exp_seeds} per production scenario')
    per_cell = args.exp_seeds * len(SCENARIOS)
    print(f'  production precision: {per_cell} simulations per cell, '
          f'completion +/-{se_pct(per_cell // 2, per_cell):.0f} points')
    if done < len(cells):
        print()
        print('  *** cells without a calibration table did not finish. The')
        print('  *** tables below pool a different mix and are not comparable.')


def rate(row, key):
    # control has no long shedders, so nucleotide_substitution_rate_long is
    # written as an empty cell. row.get(key, 'nan') does NOT cover that: the
    # key is present, its value is '', and float('') raises.
    try:
        return float(row.get(key) or 'nan')
    except (TypeError, ValueError):
        return float('nan')


def section_rates(cells, data):
    heading('1. DOES THE CALIBRATED RATE MOVE WITH THE CONFIGURATION')
    print('  The question the earlier grids assumed away. If these columns are')
    print('  flat, calibrating once was fine and #905-#908 stand. If they move,')
    print('  every configuration measured against a single frozen table was')
    print('  measuring the wrong thing.')
    print()
    print(f"{'cell':<29}{'scenario':<12}{'NSR std':>12}{'NSR long':>12}")
    for cell in cells:
        rows = data[cell]['table']
        if not rows:
            print(f'{cell_label(cell)[:27]:<29}{"-- no table --":<12}')
            continue
        for row in rows:
            long = rate(row, 'nucleotide_substitution_rate_long')
            print(f'{cell_label(cell)[:27]:<29}{row.get("scenario_name", "?")[:10]:<12}'
                  f'{rate(row, "nucleotide_substitution_rate"):>12.6f}'
                  + (f'{long:>12.6f}' if long == long else f'{"--":>12}'))
    print()
    print('  Spread per scenario, across every cell that produced a table:')
    by_scenario = defaultdict(lambda: ([], []))
    for cell in cells:
        for row in data[cell]['table'] or []:
            name = row.get('scenario_name', '?')
            # collected independently: paired appends dropped a valid std rate
            # whenever the long one was blank, which is every control row.
            for slot, key in ((0, 'nucleotide_substitution_rate'),
                              (1, 'nucleotide_substitution_rate_long')):
                value = rate(row, key)
                if value == value:
                    by_scenario[name][slot].append(value)
    print(f"{'scenario':<14}{'NSR std min':>14}{'max':>12}{'spread':>10}"
          f"{'NSR long min':>15}{'max':>12}{'spread':>10}")
    for name, (std, long) in sorted(by_scenario.items()):
        if not std:
            continue
        s_spread = (max(std) / min(std) - 1) * 100 if min(std) else float('nan')
        line = (f'{name[:12]:<14}{min(std):>14.6f}{max(std):>12.6f}'
                f'{s_spread:>9.0f}%')
        if long:
            l_spread = ((max(long) / min(long) - 1) * 100 if min(long)
                        else float('nan'))
            line += (f'{min(long):>15.6f}{max(long):>12.6f}'
                     f'{l_spread:>9.0f}%')
        else:
            line += f'{"-- no long shedders --":>37}'
        print(line)
    print()
    print('  A spread of a few percent is calibration noise at these seed')
    print('  counts. Tens of percent is the configuration moving the rate.')


def section_configuration(cells, data):
    heading('2. THE DECISION TABLE          now measured at each cell\'s OWN rate')
    print(f"{'cell':<29}{'n':>5}{'horizon':>9}{'saturated':>11}"
          f"{'phi=1 day':>11}{'window':>9}{'never phi':>11}")
    for cell in cells:
        metrics = data[cell]['metrics']
        if not metrics or not metrics['n']:
            print(f'{cell_label(cell)[:27]:<29}{"-- no production --":>5}')
            continue
        m = metrics
        phi = sorted(m['phi_days'])
        window = sorted(m['window'])
        print(f'{cell_label(cell)[:27]:<29}{m["n"]:>5}'
              f'{100.0 * m["full"] / m["n"]:>8.0f}%'
              f'{m["saturated"]:>11}'
              f'{(phi[len(phi) // 2] if phi else float("nan")):>11.0f}'
              f'{(window[len(window) // 2] if window else float("nan")):>9.0f}'
              f'{m["phi_never"]:>11}')
    print()
    print('  Same measurement as test_config_grid (its measure() is imported,')
    print('  not reimplemented), so these numbers are directly comparable with')
    print('  #905-#908 -- which is the whole point of running it this way.')


def section_cost(cells, data, outcomes, wall):
    heading('3. COST')
    if wall:
        print(f'  grid wall clock : {hms(wall)}')
    if outcomes:
        times = sorted(v[1] for v in outcomes.values())
        if times:
            print(f'  per-cell pipeline: median {hms(times[len(times) // 2])}, '
                  f'max {hms(times[-1])}')
    print()
    print('  A cell is cal_1 -> cal_2 -> production, strictly sequential; the')
    print('  cells themselves ran concurrently.')


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--populations', type=int, nargs='+',
                        default=[1000, 3000, 5000])
    parser.add_argument('--i0-fractions', type=float, nargs='+',
                        default=[0.01, 0.03])
    parser.add_argument('--r-values', type=float, nargs='+',
                        default=[1.05, 1.06, 1.07])
    parser.add_argument('--exp-base', type=int, default=1,
                        help='experiment number of the first cell; the rest '
                             'follow in order (default %(default)s). Bump it '
                             'if those numbers are already used on the '
                             'cluster -- a collision stops the stage with '
                             '"You already run an experiment with the same '
                             'name!" rather than overwriting anything.')
    parser.add_argument('--cells', type=int, nargs='+', default=None,
                        metavar='N',
                        help='run only these cells, numbered from 1 in the '
                             'order listed by --dry-run. The others keep their '
                             'experiment numbers, so a trial cell does not '
                             'renumber the grid.')
    parser.add_argument('--cal-seeds', type=int, default=10)
    parser.add_argument('--exp-seeds', type=int, default=10)
    parser.add_argument('--nsr-steps', type=int, default=5,
                        help='points in each NSR sweep (default %(default)s; '
                             'the pipeline ships 10)')
    parser.add_argument('--cal2-final-time', type=float, default=730.0,
                        help='cal_2 window (default %(default)s; the pipeline '
                             'ships 365, which at these configurations ends '
                             'before phi=1)')
    parser.add_argument('--consensus', nargs='+', default=['distribution'])
    parser.add_argument('--runner', default='slurm',
                        choices=['serial', 'multiprocessing', 'slurm'])
    parser.add_argument('--concurrency', type=int, default=6,
                        help='cells running at once (default %(default)s)')
    parser.add_argument('--cap', type=int, default=200,
                        help='total Slurm tasks released across all cells')
    parser.add_argument('--target-osr-std', type=float, default=0.0013)
    parser.add_argument('--target-osr-long', type=float, default=0.00205)
    parser.add_argument('--analyse-only', action='store_true')
    parser.add_argument('--rerun', action='store_true',
                        help='resume part-finished cells: every simulation '
                             'that already carries .completed is kept and '
                             'skipped, everything else is cleared back to '
                             'unstarted and run again, including anything '
                             'marked .failed. Without this a cell that already '
                             'started stops with "You already run an '
                             'experiment with the same name!".')
    parser.add_argument('--skip-completed', action='store_true',
                        help='skip cells that already have a calibration table '
                             'and a full production run. This is what makes a '
                             'multi-day sequential run resumable: re-launch '
                             'the same command and it picks up where it '
                             'stopped.')
    parser.add_argument('--dry-run', action='store_true',
                        help='print each cell\'s command and environment '
                             'without running anything')
    add_slurm_resource_args(parser)
    args = parser.parse_args()

    all_cells = build_cells(args.populations, args.i0_fractions, args.r_values)
    numbers = assign_numbers(all_cells, args.exp_base)
    if args.cells:
        chosen = set(args.cells)
        bad = chosen - set(range(1, len(all_cells) + 1))
        if bad:
            raise SystemExit(f'--cells out of range: {sorted(bad)} '
                             f'(1..{len(all_cells)})')
        cells = [c for i, c in enumerate(all_cells, 1) if i in chosen]
    else:
        cells = all_cells
    per_cell = (3 * args.nsr_steps * args.cal_seeds       # cal_1 groups
                + args.nsr_steps * args.cal_seeds          # cal_2 sweep
                + len(SCENARIOS) * args.exp_seeds * len(args.consensus))

    rule('=')
    print('CALIBRATED CONFIGURATION GRID')
    rule('=')
    print(f'cells        : {len(cells)}  '
          f'({len(args.populations)} populations x '
          f'{len(args.i0_fractions)} seedings x {len(args.r_values)} R)')
    print(f'per cell     : cal_1 + cal_2 + production, ~{per_cell} tasks')
    print(f'total        : ~{per_cell * len(cells)} tasks')
    print(f'R_long       : set equal to R in every cell')
    print(f'cal_2 window : {args.cal2_final_time:g} days   '
          f'NSR steps: {args.nsr_steps}')
    print(f'exp numbers  : {args.exp_base}..{args.exp_base + len(all_cells) - 1}'
          f'   one per cell, in order')
    if args.cells:
        print(f'running      : cells {sorted(args.cells)} of {len(all_cells)}  '
              f'(numbers {sorted(numbers[c] for c in cells)})')
    print()

    if args.dry_run:
        heading('DRY RUN')
        cap_per_cell = max(1, args.cap // max(1, min(args.concurrency,
                                                     len(cells))))
        for index, cell in enumerate(all_cells, 1):
            exp_num = numbers[cell]
            env = cell_env(cell, args, cap_per_cell)
            mark = ' ' if (not args.cells or index in set(args.cells)) else '-'
            print(f'{mark} {index:>3}. {cell_label(cell)}  -> #{exp_num}')
            print('     ' + ' '.join(cell_command(cell, exp_num, args)))
            print('     ' + '  '.join(
                f'{k}={env[k]}' for k in sorted(env)
                if k.startswith('SIMPLICITY_CELL')
                or k in ('SIMPLICITY_RESUME',
                         'SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM')))
        rule('=')
        return

    outcomes, wall = {}, None
    pending = cells
    if args.skip_completed:
        pending = [c for c in cells if not cell_is_complete(c, numbers[c], args)]
        done = len(cells) - len(pending)
        if done:
            print(f'[grid] {done} cell(s) already complete, skipping them')
    if not args.analyse_only and pending:
        heading('RUNNING')
        print(f'  {args.concurrency} cells at once, '
              f'{max(1, args.cap // max(1, min(args.concurrency, len(cells))))} '
              f'tasks released per cell\n')
        start = time.monotonic()
        outcomes = run_cells(pending, numbers, args)
        wall = time.monotonic() - start
        print()

    data = collect(cells, numbers, args)
    section_completeness(cells, data, outcomes, args)
    section_rates(cells, data)
    section_configuration(cells, data)
    section_cost(cells, data, outcomes, wall)
    rule('=')

    out = os.path.join('Data', f'calibrated_grid_report_{args.exp_base}.txt')
    try:
        os.makedirs('Data', exist_ok=True)
        with open(out, 'w') as handle:
            json.dump({cell_short(c): {'exp_num': numbers[c],
                                       'n_production': data[c]['n_production']}
                       for c in cells}, handle, indent=1)
        print(f'[grid] per-cell index written to {out}')
    except OSError as exc:
        print(f'[grid] could not write the index: {exc}')


if __name__ == '__main__':
    main()
