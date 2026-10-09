#!/usr/bin/env python3
'''The multiprocessing runner, which the refactor changed and nothing has run.

    python tests/test_multiprocessing_runner.py

WHY

Stage 5 rewrote all four runners to address a repeat by (group, index) instead
of a filesystem path. serial is exercised by every other test here; slurm is
covered by tests/fakeslurm and verified on the cluster. multiprocessing was
migrated at the same time and has not been run since.

It is the one runner that crosses a PROCESS boundary without Slurm: it submits
`(run_seeded_simulation, experiment_name, group, index)` to a
ProcessPoolExecutor, so the arguments are pickled and the child re-imports
everything. A repeat identified by a path survived that trivially; one
identified by (group, index) has to be resolvable in a fresh interpreter, from
files, with no inherited state. That is worth checking once.

What it asserts: every repeat ran exactly once, in its own process, and each
produced output under its own group and seed.
'''
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import simplicity.dir_manager as dm
import simplicity.jobs as jobs
import simplicity.settings_manager as sm

failures = []


def check(label, got, want):
    ok = got == want
    print(f"  [{'ok ' if ok else 'FAIL'}] {label}: {got!r}")
    if not ok:
        failures.append(f'{label}: got {got!r}, want {want!r}')


def main():
    root = tempfile.mkdtemp(prefix='mp_runner_')
    name = 'mp_probe'
    data_dir = os.path.join(root, 'Data')
    saved = os.getcwd()
    try:
        # the child processes resolve Data/ from their cwd, as dir_manager does
        os.chdir(root)
        # before anything resolves a path into it: get_reference() writes
        # reference.txt into the data dir on first import
        os.makedirs(data_dir, exist_ok=True)
        dm.set_data_dir(data_dir)
        sm._data_dir = data_dir
        os.environ['SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_MULTIPROCESS'] = '2'

        import simplicity.runme as runme
        import simplicity.runners.multiprocessing as mp

        print('running two groups through a ProcessPoolExecutor ...\n')
        groups = [{'name': 'alpha', 'R': [1.1, 1.4]},
                  {'name': 'beta', 'R': [1.7]}]
        runme.run_experiment(
            name,
            lambda: ({'_scenario_groups': groups},
                     {'population_size': 40, 'infected_individuals_at_start': 4,
                      'final_time': 12}, 2),
            simplicity_runner=mp, archive_experiment=False)

        print('\nevery repeat ran, in its own process')
        expected = jobs.all_repeats(name)
        check('six repeats were defined', len(expected), 6)

        produced = []
        for group, record in expected:
            out = os.path.join(dm.get_experiment_output_dir(name), group,
                               record['stem'], jobs.repeat_label(record))
            produced.append(os.path.isfile(os.path.join(out, 'final_time.csv')))
        check('each one wrote its output under its own group and seed',
              produced, [True] * 6)

        print('\nand the tree is partitioned by group')
        check('groups on disk', sorted(os.listdir(
            dm.get_experiment_output_dir(name))), ['alpha', 'beta'])
        check('alpha has two simulations',
              len(dm.get_simulation_output_dirs(name, group='alpha')), 2)
        check('beta has one',
              len(dm.get_simulation_output_dirs(name, group='beta')), 1)
        check('and the swept values resolve',
              sorted({sm.get_parameter_value_from_simulation_output_dir(s, 'R')
                      for s in dm.get_simulation_output_dirs(name)}),
              [1.1, 1.4, 1.7])
    finally:
        os.chdir(saved)
        shutil.rmtree(root, ignore_errors=True)

    print('\n' + '=' * 70)
    if failures:
        print(f'{len(failures)} FAILURE(S):')
        for line in failures:
            print(f'  {line}')
    else:
        print('the multiprocessing runner resolves repeats across processes')
    sys.exit(1 if failures else 0)


if __name__ == '__main__':
    main()
