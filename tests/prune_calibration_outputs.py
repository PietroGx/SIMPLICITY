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

PRODUCTION IS NEVER TOUCHED. Its figures read lineage_frequency.csv, and
PROD_EXP_NAME is a strict prefix of the unbound calibration names
("impact_long_shedders_unbound" vs "..._cal_long"), so the match is anchored on
the full name followed by "_#<number>" rather than on a prefix.

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


def calibration_prefixes():
    """Experiment-name prefixes that are calibration stages, from the configs."""
    names = set()
    try:
        import impact_long_shedders_unbound_config as unbound
        names.update(filter(None, (getattr(unbound, n, None)
                                   for n in ('LONG_NSR_EXP_NAME',
                                             'STD_NSR_EXP_NAME'))))
    except Exception as exc:
        print(f'[warn] could not read the unbound config: {exc}')
    try:
        import impact_long_shedders_config as bound
        names.update(filter(None, (getattr(bound, n, None)
                                   for n in ('LONG_NSR_EXP_NAME',
                                             'STD_NSR_EXP_NAME'))))
    except Exception as exc:
        print(f'[warn] could not read the bound config: {exc}')
    # the bound pipeline's standard sweep predates that constant
    names.add('impact_long_shedders_calibration_std_nsr')
    return sorted(names)


def matching_experiments(data_dir, prefixes):
    """Directories named <calibration prefix>_#<number>, and nothing else."""
    patterns = [re.compile(rf'^{re.escape(p)}_#\d+$') for p in prefixes]
    found = []
    if not os.path.isdir(data_dir):
        return found
    for entry in sorted(os.listdir(data_dir)):
        path = os.path.join(data_dir, entry)
        if os.path.isdir(path) and any(p.match(entry) for p in patterns):
            found.append(path)
    return found


def prunable_files(experiment_dir):
    out = []
    for dirpath, _, files in os.walk(os.path.join(experiment_dir, '04_Output')):
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

    prefixes = calibration_prefixes()
    print('calibration experiment names (from the pipeline config):')
    for prefix in prefixes:
        print(f'    {prefix}_#<n>')

    experiments = matching_experiments(args.data_dir, prefixes)
    if not experiments:
        print(f'\nNo calibration experiments under {args.data_dir}/.')
        return

    print(f'\nmatched {len(experiments)} calibration experiment(s); '
          f'everything else, production included, is left alone\n')
    print(f"{'experiment':<52}{'files':>8}{'size':>12}")
    total_files = total_bytes = 0
    per_experiment = []
    for path in experiments:
        files = prunable_files(path)
        size = sum(s for _, s in files)
        per_experiment.append((path, files))
        total_files += len(files)
        total_bytes += size
        print(f'{os.path.basename(path)[:50]:<52}{len(files):>8}'
              f'{human(size):>12}')

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
