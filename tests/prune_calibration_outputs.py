#!/usr/bin/env python3
'''Delete the outputs calibration never reads, from calibration runs already on
disk.

v2.4.65 stopped writing them; this applies the same decision retroactively.

    python tests/prune_calibration_outputs.py              # report only
    python tests/prune_calibration_outputs.py --apply      # delete

Dry run by default, and it prints every experiment it matched before touching
anything -- read that list first. Deletion is irreversible and this operates on
the real Data/, not a copy.

WHAT IT DELETES, AND WHY THAT IS SAFE

Only inside CALIBRATION experiments, whose names come from the pipeline config
rather than from a pattern guessed here, and only these three files:

    lineage_frequency.csv       nothing in the calibration or sanity path opens it
    simulation_trajectory.csv   likewise
    sequencing_data.csv         no reader anywhere in the repo -- sequencing.py
                                reconstructs those rows from individuals_data
                                plus phylogenetic_data

What stays, because the calibration path does read it: final_time.csv and
sequencing_data_regression.csv for both OSR fits, individuals_data.csv and
phylogenetic_data.csv for the intra-host clock through
evolutionary_rate.extract_ih_regression_data.

THE EXPERIMENT PHASE IS NEVER TOUCHED. Its figures read
lineage_frequency.csv, and EXP_NAME is a strict prefix of the calibration
names ("impact_long_shedders" vs "..._cal_long"), so the match is anchored on
the full name followed by "_#<number>" rather than on a prefix. That prefix
relationship is why it is anchored, and it got tighter in v2.4.91 when the
"_unbound" infix went away.

Safe to run on a part-finished experiment too: a resumed run re-runs whatever
lacks a .completed signal, and these files are not read on the way there.
'''
import argparse
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'scripts', 'experiments'))

PRUNABLE = ('lineage_frequency.csv', 'simulation_trajectory.csv',
            'sequencing_data.csv')


# Experiment names earlier versions gave a calibration stage of its own.
# The code that wrote them is gone -- the bound pipeline in v2.4.90, the
# per-stage experiments in v2.4.94 -- but their runs are still on disk and
# still prunable, so the names stay here.
_LEGACY_CALIBRATION_PREFIXES = (
    'impact_long_shedders_calibration_lng_nsr',
    'impact_long_shedders_calibration_std_nsr',
    'impact_long_shedders_cal_long',
    'impact_long_shedders_cal_std',
    'impact_long_shedders_unbound_cal_long',
    'impact_long_shedders_unbound_cal_std',
)


def calibration_prefixes():
    """Experiment-name prefixes that were whole calibration stages."""
    return sorted(_LEGACY_CALIBRATION_PREFIXES)


def calibration_groups(experiment_dir):
    """The calibration GROUPS of a merged experiment, or [] if it has none.

    Since 8e a calibration is a set of groups inside the run's one
    experiment, not an experiment of its own, so matching on the directory
    name finds nothing. The experiment phase's groups live in the same
    directory and must NOT be pruned -- the figures read
    lineage_frequency.csv -- so this returns only the calib_* ones and the
    walk below is scoped to them.
    """
    try:
        import impact_long_shedders_config as cfg
    except Exception as exc:                                    # noqa: BLE001
        print(f'[warn] could not read the config: {exc}')
        return []
    output = os.path.join(experiment_dir, '04_Output')
    if not os.path.isdir(output):
        return []
    return [name for name in sorted(os.listdir(output))
            if os.path.isdir(os.path.join(output, name))
            and cfg.phase_of(name) in ('calib_long', 'calib_std')]


def matching_targets(data_dir):
    """Every (experiment dir, groups) worth pruning.

    groups is None for a legacy stand-alone calibration experiment -- all of
    it is calibration -- and a list for a merged experiment, where only the
    calib_* groups are.
    """
    patterns = [re.compile(rf'^{re.escape(p)}_#\d+$')
                for p in calibration_prefixes()]
    found = []
    if not os.path.isdir(data_dir):
        return found
    for entry in sorted(os.listdir(data_dir)):
        path = os.path.join(data_dir, entry)
        if not os.path.isdir(path):
            continue
        if any(p.match(entry) for p in patterns):
            found.append((path, None))
            continue
        groups = calibration_groups(path)
        if groups:
            found.append((path, groups))
    return found


def prunable_files(experiment_dir, groups=None):
    out = []
    roots = ([os.path.join(experiment_dir, '04_Output', g) for g in groups]
             if groups else [os.path.join(experiment_dir, '04_Output')])
    for root in roots:
      for dirpath, _, files in os.walk(root):
        for name in files:
            if name in PRUNABLE:
                path = os.path.join(dirpath, name)
                try:
                    out.append((path, os.path.getsize(path)))
                except OSError:
                    pass
    return out


def human(n):
    for unit in ('B', 'KB', 'MB', 'GB', 'TB'):
        if n < 1024 or unit == 'TB':
            return f'{n:.1f} {unit}'
        n /= 1024


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--data-dir', default='Data')
    parser.add_argument('--apply', action='store_true',
                        help='actually delete. Without it, nothing is removed.')
    args = parser.parse_args()

    print('a calibration is pruned either way:')
    print('  as a whole experiment, for the names older versions gave one:')
    for prefix in calibration_prefixes():
        print(f'    {prefix}_#<n>')
    print('  or as the calib_long__* / calib_std__* GROUPS of a merged '
          'experiment.')

    targets = matching_targets(args.data_dir)
    if not targets:
        print(f'\nNothing to prune under {args.data_dir}/.')
        return

    print(f'\nmatched {len(targets)} experiment(s); every other group, the '
          f'experiment phase included, is left alone\n')
    print(f"{'experiment / group':<52}{'files':>8}{'size':>12}")
    total_files = total_bytes = 0
    per_experiment = []
    for path, groups in targets:
        files = prunable_files(path, groups)
        size = sum(s for _, s in files)
        per_experiment.append((path, files))
        total_files += len(files)
        total_bytes += size
        label = os.path.basename(path)
        if groups:
            label += f'  [{len(groups)} calib group(s)]'
        print(f'{label[:50]:<52}{len(files):>8}{human(size):>12}')

    print(f'\n{"TOTAL":<52}{total_files:>8}{human(total_bytes):>12}')

    if not args.apply:
        print('\nDry run. Nothing was deleted.')
        print('Read the list above, then re-run with --apply.')
        return

    print('\ndeleting...')
    removed = freed = 0
    for path, files in per_experiment:
        for file_path, size in files:
            try:
                os.remove(file_path)
                removed += 1
                freed += size
            except OSError as exc:
                print(f'  [warn] {file_path}: {exc}')
    print(f'removed {removed} file(s), freed {human(freed)}')


if __name__ == '__main__':
    main()
