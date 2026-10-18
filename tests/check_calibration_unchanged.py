#!/usr/bin/env python3
'''Re-fit a finished calibration with the refactored code and compare it to
what that run actually produced. Run this ON THE CLUSTER, where the output is.

    python tests/check_calibration_unchanged.py --exp-num 2
    python tests/check_calibration_unchanged.py --exp-num 2 --keep
    python tests/check_calibration_unchanged.py --exp-num 2 --only 1

WHY

Stage 7c changed HOW a calibration decides which simulations belong to which
group: by the group each one was assigned to, instead of by reading
(tau_3_long, R_long) or (tau_3_long, long_shedders_ratio) back off every output
directory and bucketing on the float values.

tests/test_group_partition.py already shows the two agree on the real grids --
but it shows it on grids it built itself, with a synthetic measurement. This
shows it on an actual finished run, with the real measurement, by recomputing
the calibration and diffing against the numbers that run froze. If they match,
the refactor changed nothing that was ever right and no recalibration is owed.
If they do not, that difference is the whole finding and the frozen table it
came from is suspect.

This is the check that could not be done when 7c was deferred: there was no
Data/ to verify against.

WHAT IT TOUCHES

Recomputing overwrites the experiment's Fit_results and its calibration plot,
so both are copied aside first and PUT BACK before this exits. --keep leaves
the recomputed ones in place instead. Nothing is simulated and nothing is
submitted; it only re-reads output that is already on disk.

The standard-NSR target is not recorded in the frozen table, so --target-osr-std
defaults to the pipeline's own default and is printed. If that run used another
value, pass it, or stage 2's comparison is against the wrong target.
'''
import argparse
import os
import shutil
import sys
import tempfile

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'scripts'))
sys.path.insert(0, os.path.join(REPO, 'scripts', 'experiments'))

import simplicity.dir_manager as dm
from impact_long_shedders_config import (
    EXP_NAME, calib_std_group, table_path as frozen_table_path)

failures = []
skipped = []

# One experiment since 8e; the two calibrations are GROUPS inside it. Where
# the frozen table lives is NOT spelled here -- the config's table_path owns
# that, and since v2.4.111 it nests inside the run's own directory.
NAMES = {
    'cal_1': EXP_NAME,
    'cal_2': EXP_NAME,
}
TABLE_FILENAME = 'nsr_calibration_table.csv'


def compare(label, got, want, rtol):
    if want == 0:
        ok = got == 0
        delta = 'n/a'
    else:
        delta = f'{100 * (got - want) / abs(want):+.4f}%'
        ok = abs(got - want) / abs(want) <= rtol
    print(f"  [{'ok ' if ok else 'FAIL'}] {label:<22} "
          f"was {want:.10g}   now {got:.10g}   ({delta})")
    if not ok:
        failures.append(f'{label}: was {want!r}, recomputed {got!r} ({delta})')


class Preserved:
    '''Copy the directories a re-fit writes into, and put them back.'''

    def __init__(self, *directories):
        self.pairs = []
        self.backup = tempfile.mkdtemp(prefix='calib_backup_')
        for index, directory in enumerate(directories):
            if directory and os.path.isdir(directory):
                saved = os.path.join(self.backup, str(index))
                shutil.copytree(directory, saved)
                self.pairs.append((directory, saved))

    def restore(self):
        for directory, saved in self.pairs:
            shutil.rmtree(directory, ignore_errors=True)
            shutil.copytree(saved, directory)
        self.discard()

    def discard(self):
        shutil.rmtree(self.backup, ignore_errors=True)


def experiment_exists(name):
    try:
        return bool(dm.get_groups(name))
    except Exception:                                           # noqa: BLE001
        return False


# ===========================================================================
def stage_1(names, exp_num, rtol, keep):
    '''cal_1: the calibrated NSR_long per long-shedder cohort.'''
    import long_nsr_calibration_plot as plot

    name = f"{names['cal_1']}_#{exp_num}"
    print(f'\n=== stage 1: {name}')
    if not experiment_exists(name):
        skipped.append(f'stage 1: no experiment {name}')
        print('  SKIPPED -- no such experiment on this machine')
        return

    frozen_path = plot._calibrated_long_nsr_path(name)
    if not os.path.isfile(frozen_path):
        skipped.append(f'stage 1: {name} never wrote a calibrated_long_nsr.csv')
        print(f'  SKIPPED -- {frozen_path} does not exist, so there is '
              'nothing to compare against')
        return

    before = pd.read_csv(frozen_path)
    target = float(before['target_osr_long'].iloc[0])
    model = str(before['model_type'].iloc[0])
    print(f'  recorded target_osr_long={target}  model={model}  '
          f'{len(before)} cohort(s)')
    print(f'  groups on disk: {dm.get_groups(name)}')

    preserved = Preserved(dm.get_experiment_fit_result_dir(name),
                          dm.get_experiment_plots_dir(name))
    try:
        results = plot.plot_and_fit_long_nsr_calibration(
            name, target, model_type=model)
    except Exception as exc:                                    # noqa: BLE001
        failures.append(f'stage 1 re-fit raised {type(exc).__name__}: {exc}')
        print(f'  [FAIL] the re-fit raised {type(exc).__name__}: {exc}')
        preserved.restore()
        return

    print()
    old = {(round(float(r['tau_3_long']), 3), round(float(r['R_long']), 3)):
           float(r['nucleotide_substitution_rate_long'])
           for _, r in before.iterrows()}
    if sorted(old) != sorted(results):
        failures.append(f'stage 1: cohorts changed -- was {sorted(old)}, '
                        f'now {sorted(results)}')
        print('  [FAIL] the set of cohorts changed')
        print(f'         was {sorted(old)}')
        print(f'         now {sorted(results)}')
    for key in sorted(set(old) & set(results)):
        compare(f'tau={key[0]} R={key[1]}', results[key], old[key], rtol)

    if keep:
        preserved.discard()
        print('  (recomputed Fit_results and plot KEPT)')
    else:
        preserved.restore()
        print('  (Fit_results and plot restored to what the run produced)')


# ===========================================================================
def stage_2(names, exp_num, target_osr_std, rtol, keep):
    '''cal_2: the calibrated standard NSR, against the frozen table.'''
    name = f"{names['cal_2']}_#{exp_num}"
    table_path = frozen_table_path(exp_num)
    print(f'\n=== stage 2: {name}')
    if not experiment_exists(name):
        skipped.append(f'stage 2: no experiment {name}')
        print('  SKIPPED -- no such experiment on this machine')
        return
    if not os.path.isfile(table_path):
        skipped.append(f'stage 2: no frozen table at {table_path}')
        print(f'  SKIPPED -- no frozen table at {table_path}, so there is '
              'nothing to compare against')
        return

    table = pd.read_csv(table_path)
    model = str(table['model_type'].iloc[0]) if 'model_type' in table else 'exp'
    print(f'  target_osr_std={target_osr_std} (NOT recorded in the table -- '
          'pass --target-osr-std if this run used another)')
    print(f'  model={model}   groups on disk: {dm.get_groups(name)}')

    preserved = Preserved(dm.get_experiment_fit_result_dir(name),
                          dm.get_experiment_plots_dir(name))
    try:
        import impact_long_shedders_cal_2 as cal_2
        single = cal_2.compute_standard_nsr(name, target_osr_std,
                                            model_type=model,
                                            groups=[calib_std_group()])
        # one standard rate for every scenario, by construction: the same
        # number is printed once per row below
        recomputed = {row['scenario_name']: single
                      for _, row in table.iterrows()}
    except Exception as exc:                                    # noqa: BLE001
        failures.append(f'stage 2 re-fit raised {type(exc).__name__}: {exc}')
        print(f'  [FAIL] the re-fit raised {type(exc).__name__}: {exc}')
        preserved.restore()
        return

    print()
    for _, row in table.iterrows():
        scenario = row['scenario_name']
        if scenario not in recomputed:
            failures.append(f'stage 2: {scenario} was not recomputed at all')
            print(f'  [FAIL] {scenario}: not recomputed')
            continue
        compare(scenario, recomputed[scenario],
                float(row['nucleotide_substitution_rate']), rtol)

    if keep:
        preserved.discard()
        print('  (recomputed Fit_results and plot KEPT)')
    else:
        preserved.restore()
        print('  (Fit_results and plot restored to what the run produced)')


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--exp-num', type=int, required=True)
    parser.add_argument('--target-osr-std', type=float, default=0.0013,
                        help='the pipeline default; not recorded in the table')
    parser.add_argument('--rtol', type=float, default=1e-9,
                        help='relative tolerance. The default is strict on '
                             'purpose: a pure repartition that still reads the '
                             'same rows should reproduce the fit bit for bit, '
                             'so anything above float noise is a real change.')
    parser.add_argument('--keep', action='store_true',
                        help='leave the recomputed Fit_results and plots in '
                             'place instead of restoring the originals')
    parser.add_argument('--only', choices=['1', '2'], default=None)
    args = parser.parse_args()

    if not os.path.isdir('Data'):
        print('error: no Data/ here. Run this from the repo root, on the '
              'machine that holds the run.', file=sys.stderr)
        sys.exit(2)

    print(f'run #{args.exp_num}, rtol={args.rtol:g}')

    if args.only != '2':
        stage_1(NAMES, args.exp_num, args.rtol, args.keep)
    if args.only != '1':
        stage_2(NAMES, args.exp_num, args.target_osr_std, args.rtol, args.keep)

    print('\n' + '=' * 70)
    for line in skipped:
        print(f'  SKIPPED  {line}')
    if failures:
        print(f'{len(failures)} DIFFERENCE(S) -- the refactor changed a result:')
        for line in failures:
            print(f'  {line}')
    elif skipped and len(skipped) >= 2:
        print('nothing was compared: neither stage had both a finished run and '
              'the\nnumbers it froze. This proves nothing -- point it at a run '
              'that completed.')
        sys.exit(1)
    else:
        print('the refactored code reproduces this run\'s calibration exactly.\n'
              'No recalibration is owed.')
    sys.exit(1 if failures else 0)


if __name__ == '__main__':
    main()
