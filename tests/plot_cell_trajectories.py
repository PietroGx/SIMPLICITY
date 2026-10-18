#!/usr/bin/env python3
'''Infected-vs-time for every simulation of one or more configuration cells.

    cd ~/SIMPLICITY && git pull
    python tests/plot_cell_trajectories.py --exp-nums 122 134 --slurm

The decision table gives a cell one number per column, and two cells with the
same horizon fraction can look nothing alike: an epidemic that runs flat to
day 1095 and one that spikes at day 300 and crawls both count as "reached the
horizon". This draws the curves the table summarises.

Reads Data/<exp>/04_Output/<group>/<simulation>/<seed>/simulation_trajectory.csv
-- the same file test_config_grid.measure() reads, and phi=1 is computed here
by the same rule ((diagnosed+recovered)/N >= 1), so the picture and the table
cannot drift apart.

N and I0 are read off the trajectory (susceptibles+infected at t=0), not from
test_calibrated_grid's axis defaults: those defaults are 3x2x3 and #101-#136
were run as 3x3x4, so they would mislabel every cell.

Two figures per cell, into Data/grid_trajectories/:
  *_by_scenario.png       5 panels, one per scenario, 10 seeds overlaid
  *_per_simulation.png    one small panel per simulation, 5 x 10
'''
import argparse
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
SCRIPTS = os.path.join(REPO, 'scripts', 'experiments')
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)
sys.path.insert(0, SCRIPTS)

import impact_long_shedders_config as cfg
from impact_long_shedders_config import SCENARIOS, DEFAULT_CONSENSUS

HORIZON = 1095.0
OUT_DIR = os.path.join('Data', 'grid_trajectories')

# extrande.check_stop_conditions has FIVE conditions; test_config_grid.measure()
# classifies four and calls the rest "other", so the decision table cannot say
# which of the last two ended a run. Both are computational, not epidemic:
#
#   len(population.reservoir_i) < 1000   extrande.py:257
#   time.time() - start > max_runtime    extrande.py:259
#
# The reservoir is a fixed pool of 100000 individuals, NOT scaled by
# population_size (population.py:78, never overridden -- population.py:704).
# It starts at 100000 - N and loses exactly one individual per diagnosis
# (population_model.py:65 and :123, same function, no return between them) and
# one per recovery (population.py:455 and :460). So the number drawn is
# diagnosed + recovered, both of which the trajectory carries, and the stop
# fires at (diagnosed + recovered) = 100000 - N - 1000:
#
#   N=5000 -> 94000 events -> phi = 18.8      N=3000 -> 96000 -> phi = 32.0
#
# Five conditions, four ruled out, so the fifth is named rather than guessed.
RESERVOIR = 100000
RESERVOIR_FLOOR = 1000

# The experiment names a cell may be on disk under: 'ils_grid' from v2.4.111,
# 'impact_long_shedders' for #101-#136, which ran before the rename.
NAME_CANDIDATES = ('ils_grid', 'impact_long_shedders')


def burnin_cutoff():
    """fig3_preprocess_data.BURNIN_CUTOFF_DAYS, read rather than re-spelled.

    Imported it would drag in _scenarios and long_shedders_preprocess through
    two sibling-directory sys.path tricks; the constant is what is wanted."""
    path = os.path.join(REPO, 'scripts', 'nature_plots', 'fig3_preprocess_data.py')
    try:
        m = re.search(r'^BURNIN_CUTOFF_DAYS\s*=\s*(\d+)', open(path).read(), re.M)
        return float(m.group(1))
    except Exception:
        return None


# ------------------------------------------------------------------ locating

def resolve_run(exp_num, consensus):
    """(experiment name, run directory) for a cell, or (None, None).

    cfg.EXP_NAME is what every path helper derives from, so rebinding it is
    how the whole directory moves -- the same thing test_calibrated_grid does
    to put the cells under ils_grid."""
    for base in NAME_CANDIDATES:
        cfg.EXP_NAME = base
        run_dir = cfg.run_dir(exp_num, consensus)
        if os.path.isdir(run_dir):
            return base, run_dir
    return None, None


def trajectories(run_dir):
    """{group: {simulation: {seed: path}}}, sorted. Calibration groups never
    appear: CALIBRATION_SKIPPED_OUTPUTS switches the trajectory off for them,
    so only the experiment phase writes one."""
    found = {}
    root = os.path.join(run_dir, '04_Output')
    for dirpath, _, files in os.walk(root):
        if 'simulation_trajectory.csv' not in files:
            continue
        seed = os.path.basename(dirpath)
        sim = os.path.basename(os.path.dirname(dirpath))
        group = os.path.basename(os.path.dirname(os.path.dirname(dirpath)))
        found.setdefault(group, {}).setdefault(sim, {})[seed] = \
            os.path.join(dirpath, 'simulation_trajectory.csv')
    return found


def group_order(groups):
    """Scenario order as the config declares it, with anything unrecognised
    after it rather than dropped."""
    rank = {s['name']: i for i, s in enumerate(s for s in SCENARIOS)}
    def key(g):
        tail = g.split('__')[-1]
        return (rank.get(tail, len(rank)), g)
    return sorted(groups, key=key)


def cell_params(run_dir):
    """{'R', 'population_size'} from the experiment record, missing keys left
    out. Read straight off the JSON: the settings_manager accessors resolve
    paths through dir_manager's absolute _data_dir, which is not where a run
    uploaded from the cluster sits.

    population_size matters more than the label. phi is (diagnosed+recovered)/N
    and nothing in the trajectory states N, so taking it from the record rather
    than inferring it is what keeps the dotted line honest -- the same reason
    test_config_grid.measure() demands pop_size as an argument."""
    path = os.path.join(run_dir, '01_Settings', 'settings.json')
    out = {}
    try:
        sims = json.load(open(path))['simulations']
    except Exception:
        return out
    for key in ('R', 'population_size'):
        values = {s['parameters'].get(key) for s in sims}
        values.discard(None)
        if len(values) == 1:
            out[key] = values.pop()
    return out


# ------------------------------------------------------------------ reading

def daily(path, pop_size=None):
    """One run, thinned to a value per day: (days, infected, record).

    The trajectory is one row per event -- up to ~1M rows and 47 MB at
    N=5000 -- and a 1095-point step function draws the same curve. Each day
    takes the last row at or before it, so nothing is interpolated across a
    gap and the series stops where the run stopped."""
    import numpy as np
    import pandas as pd
    df = pd.read_csv(path)
    t = df['time'].to_numpy(dtype=float)
    if t.size == 0:
        return None
    # t=0 has nobody diagnosed or recovered yet, so the two compartments
    # account for everyone -- but the record is believed over the inference.
    seen = float(df['susceptibles'].iloc[0] + df['infected'].iloc[0])
    pop = float(pop_size) if pop_size else seen
    phi = (df['diagnosed'].to_numpy(dtype=float)
           + df['recovered'].to_numpy(dtype=float)) / pop if pop else None
    hit = np.flatnonzero(phi >= 1.0) if phi is not None else np.array([], int)

    days = np.arange(0.0, np.floor(t[-1]) + 1.0)
    idx = np.clip(np.searchsorted(t, days, side='right') - 1, 0, t.size - 1)
    infected = df['infected'].to_numpy(dtype=float)[idx]

    last = df.iloc[-1]
    drawn = float(last['diagnosed']) + float(last['recovered'])
    left = RESERVOIR - pop - drawn
    record = {
        'reservoir_left': left,
        'drawn': drawn,
        'pop': pop,
        'pop_seen': seen,
        'i0': float(df['infected'].iloc[0]),
        'end': float(t[-1]),
        'peak': float(df['infected'].max()),
        'peak_day': float(t[int(np.argmax(df['infected'].to_numpy()))]),
        'phi_day': float(t[hit[0]]) if hit.size else None,
        'full': float(t[-1]) >= HORIZON,
        # extrande.check_stop_conditions, in its order and complete
        'stop': ('horizon' if float(t[-1]) >= HORIZON else
                 'died' if float(last['infected']) == 0 else
                 'saturated' if float(last['susceptibles']) == 0 else
                 'no infectious' if float(last['infectious']) == 0
                 and float(t[-1]) > 60.0 else
                 'reservoir' if left < RESERVOIR_FLOOR else 'max runtime'),
    }
    return days, infected, record


def load_cell(paths_by_group, pop_size=None):
    """{group: [(seed label, days, infected, record), ...]}"""
    out = {}
    for group in group_order(paths_by_group):
        runs = []
        for sim in sorted(paths_by_group[group]):
            for seed in sorted(paths_by_group[group][sim]):
                got = daily(paths_by_group[group][sim][seed], pop_size)
                if got is None:
                    continue
                runs.append((seed, got[0], got[1], got[2]))
                print(f'    {group:<20}{seed:<12}'
                      f'end {got[2]["end"]:>7.0f}  peak {got[2]["peak"]:>6.0f}'
                      f'  drawn {got[2]["drawn"]:>7.0f}'
                      f'  reservoir left {got[2]["reservoir_left"]:>7.0f}'
                      f'  {got[2]["stop"]}')
        out[group] = runs
    return out


# ----------------------------------------------------------------- drawing

def short(group):
    return group.split('__')[-1]


def draw_by_scenario(cell, title, out_path, cutoff):
    import numpy as np
    import matplotlib.pyplot as plt
    groups = list(cell)
    fig, axes = plt.subplots(len(groups), 1, figsize=(11, 2.1 * len(groups)),
                             sharex=True, sharey=True)
    axes = np.atleast_1d(axes)
    for ax, group in zip(axes, groups):
        runs = cell[group]
        for _, days, infected, rec in runs:
            ax.plot(days, infected, lw=0.8, alpha=0.75)
            ax.plot(rec['end'], infected[-1], 'k.', ms=3)
        if cutoff is not None:
            ax.axvline(cutoff, color='k', ls='--', lw=0.8)
        phis = [r['phi_day'] for _, _, _, r in runs if r['phi_day'] is not None]
        if phis:
            ax.axvline(float(np.median(phis)), color='tab:red', ls=':', lw=1.0)
        ends = [r['end'] for _, _, _, r in runs]
        full = sum(1 for _, _, _, r in runs if r['full'])
        ax.set_ylabel('infected')
        ax.text(0.995, 0.93,
                f'{short(group)}   n={len(runs)}   horizon {full}/{len(runs)}'
                f'   median end {np.median(ends):.0f}d'
                + (f'   median phi=1 {np.median(phis):.0f}d' if phis
                   else '   phi=1 never'),
                transform=ax.transAxes, ha='right', va='top', fontsize=8)
    axes[-1].set_xlabel('time (days)')
    axes[0].set_xlim(0, HORIZON)
    axes[0].set_ylim(bottom=0)
    fig.suptitle(title + '   (dashed: burn-in cutoff, dotted: median phi=1, '
                 'dot: run end)', fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    print(f'  wrote {out_path}')


def draw_per_simulation(cell, title, out_path, cutoff):
    import numpy as np
    import matplotlib.pyplot as plt
    groups = list(cell)
    cols = max(len(cell[g]) for g in groups)
    fig, axes = plt.subplots(len(groups), cols, sharex=True, sharey=True,
                             figsize=(1.35 * cols, 1.5 * len(groups)))
    axes = np.atleast_2d(axes)
    for row, group in enumerate(groups):
        for col in range(cols):
            ax = axes[row][col]
            ax.set_xticks([0, 500, 1000])
            ax.tick_params(labelsize=6)
            if col >= len(cell[group]):
                ax.axis('off')
                continue
            seed, days, infected, rec = cell[group][col]
            ax.plot(days, infected, lw=0.8, color='tab:blue')
            if cutoff is not None:
                ax.axvline(cutoff, color='k', ls='--', lw=0.6)
            ax.set_title(f'{rec["end"]:.0f}d  {rec["stop"]}', fontsize=6)
            if col == 0:
                ax.set_ylabel(short(group), fontsize=7)
    axes[0][0].set_xlim(0, HORIZON)
    axes[0][0].set_ylim(bottom=0)
    fig.suptitle(title + '   (panel title: end day and stop reason)',
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    print(f'  wrote {out_path}')


def report_phi(cell, cutoff):
    """The spread of phi=1 across the cell, not just its median.

    The decision table reports the MEDIAN phi day, which is the wrong number
    to set a burn-in cutoff against: at the median, half the runs have not
    reached phi=1 yet and their pre-selection phase is still inside the
    window. A cutoff has to clear the tail, so the quantiles are what decide
    where it goes."""
    import numpy as np
    phis = sorted(r['phi_day'] for g in cell for *_, r in cell[g]
                  if r['phi_day'] is not None)
    never = sum(1 for g in cell for *_, r in cell[g] if r['phi_day'] is None)
    if not phis:
        print(f'  phi=1: never reached by any of {never} runs')
        return
    q = np.percentile(phis, [50, 75, 90, 100])
    print(f'  phi=1 day: median {q[0]:.0f}  p75 {q[1]:.0f}  p90 {q[2]:.0f}  '
          f'max {q[3]:.0f}  (never: {never})')
    for day in (200, 250, 300, 350):
        late = sum(1 for v in phis if v > day) + never
        flag = '  <- figures\' cutoff' if cutoff and day == cutoff else ''
        print(f'    cutoff {day}: {late:>2} of {len(phis) + never} runs still '
              f'pre-phi=1 there{flag}')


# -------------------------------------------------------------------- main

def submit(args):
    os.makedirs(os.path.join('slurm_logs', 'trajectories'), exist_ok=True)
    inner = [sys.executable, os.path.abspath(__file__),
             '--exp-nums', *[str(n) for n in args.exp_nums],
             '--consensus', args.consensus]
    cmd = ['sbatch', '--job-name=traj_plot', '--cpus-per-task=1',
           '--mem=8G', '--time=01:00:00',
           '--output=slurm_logs/trajectories/traj_%j.out',
           '--error=slurm_logs/trajectories/traj_%j.err',
           '--wrap=' + ' '.join(inner)]
    print(' '.join(cmd))
    subprocess.run(cmd, check=True)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--exp-nums', type=int, nargs='+', required=True)
    parser.add_argument('--consensus', default=DEFAULT_CONSENSUS)
    parser.add_argument('--out-dir', default=OUT_DIR)
    parser.add_argument('--slurm', action='store_true',
                        help='submit this script as a Slurm job instead of '
                             'reading ~2 GB of CSV on the login node')
    args = parser.parse_args()

    if args.slurm:
        submit(args)
        return 0

    import matplotlib
    matplotlib.use('Agg')

    cutoff = burnin_cutoff()
    os.makedirs(args.out_dir, exist_ok=True)
    missing = []
    for exp_num in args.exp_nums:
        base, run_dir = resolve_run(exp_num, args.consensus)
        if run_dir is None:
            missing.append(exp_num)
            print(f'#{exp_num}: not on disk under any of {NAME_CANDIDATES}')
            continue
        print(f'#{exp_num}: {run_dir}')
        paths = trajectories(run_dir)
        if not paths:
            missing.append(exp_num)
            print(f'  no simulation_trajectory.csv under {run_dir}/04_Output')
            continue
        params = cell_params(run_dir)
        cell = load_cell(paths, params.get('population_size'))
        any_run = next((r for g in cell for *_, r in cell[g]), None)
        if any_run['pop'] != any_run['pop_seen']:
            print(f'  WARNING: settings say N={any_run["pop"]:.0f}, the '
                  f'trajectory starts with {any_run["pop_seen"]:.0f} in '
                  f'susceptibles+infected. phi is scaled by the former.')
        title = (f'{base} #{exp_num}   N={any_run["pop"]:.0f}  '
                 f'I0={any_run["i0"]:.0f}'
                 + (f'  R={params["R"]}' if 'R' in params else ''))
        tally = {}
        for g in cell:
            for *_, rec in cell[g]:
                tally[rec['stop']] = tally.get(rec['stop'], 0) + 1
        print('  stop reasons: '
              + '  '.join(f'{k} {v}' for k, v in sorted(tally.items())))
        report_phi(cell, cutoff)
        stem = os.path.join(args.out_dir, f'{base}_{exp_num}_infected')
        draw_by_scenario(cell, title, stem + '_by_scenario.png', cutoff)
        draw_per_simulation(cell, title, stem + '_per_simulation.png', cutoff)

    print(f'\nfigures in {args.out_dir}/ -- collect them with:')
    print(f'  tar czf grid_trajectories.tgz {args.out_dir}')
    return 1 if missing else 0


if __name__ == '__main__':
    sys.exit(main())
