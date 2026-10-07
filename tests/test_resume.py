#!/usr/bin/env python3
'''Resuming an experiment must keep what finished and re-run only the rest.

Without a resume, re-submitting an experiment is refused twice over:
raise_on_stale_signals stops it before sbatch, and setup_output_directory
raises "You already run an experiment with the same name!" per simulation.
Both are right by default -- which of two runs matters is not something the
code can know. SIMPLICITY_RESUME=1 (the --rerun flag) says the earlier one
matters and the gap should be filled.

What that has to get right, and what this checks:

  completed kept      a simulation carrying .completed is left alone and its
                      output is never recomputed. Cell #17 of the calibrated
                      grid stalled with 61 of 150 finished; throwing those away
                      was the first thing I built and it was wrong.

  failed retried      a .failed task is a job still to do on a resume, not a
                      settled answer, so its signals are cleared too.

  partial cleared     a simulation with .started and no .completed has a
                      half-written output directory. It is removed before the
                      re-run: leaving it mixes the abandoned attempt's files
                      with the new ones and nothing downstream can tell which
                      is which.

  refuses by default  with SIMPLICITY_RESUME unset, both guards still fire.
                      A resume must be asked for, never inferred.

    python tests/test_resume.py
'''
import os
import pathlib
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import simplicity.runners.slurm as slurm
import simplicity.output_manager as om

failures = []


def check(label, got, want):
    ok = got == want
    print(f"  [{'ok ' if ok else 'FAIL'}] {label}: {got!r}")
    if not ok:
        failures.append(f'{label}: got {got!r}, want {want!r}')


def signals(path):
    return {s for s in ('.submitted', '.released', '.started', '.completed',
                        '.failed', '.progress')
            if pathlib.Path(path + s).exists()}


def main():
    root = tempfile.mkdtemp(prefix='resume_test_')
    original = slurm.sm.get_seeded_simulation_parameters_paths
    try:
        # four simulations, one per state a stalled experiment leaves behind
        names = ['done', 'failed', 'started_only', 'never_ran']
        paths = [os.path.join(root, f'{n}.json') for n in names]
        for path in paths:
            pathlib.Path(path).write_text('{}')
            pathlib.Path(path + '.submitted').touch()
            pathlib.Path(path + '.released').touch()
        for suffix in ('.started', '.completed'):
            pathlib.Path(paths[0] + suffix).touch()
        for suffix in ('.started', '.failed', '.progress'):
            pathlib.Path(paths[1] + suffix).touch()
        pathlib.Path(paths[2] + '.started').touch()
        pathlib.Path(paths[2] + '.progress').touch()

        slurm.sm.get_seeded_simulation_parameters_paths = lambda name: paths

        print('\nprepare_resume')
        kept, to_run = slurm.prepare_resume('irrelevant')
        check('completed simulations kept', kept, 1)
        check('everything else queued to run', to_run, 3)
        check('the finished one is untouched', signals(paths[0]),
              {'.submitted', '.released', '.started', '.completed'})
        check('the failed one is cleared back to unstarted',
              signals(paths[1]), set())
        check('the half-started one is cleared', signals(paths[2]), set())
        check('the untouched one is cleared', signals(paths[3]), set())

        print('\nslurm.job skips what finished')
        os.environ['SIMPLICITY_RESUME'] = '1'
        check('resuming() reads the environment', slurm.resuming(), True)
        check('a completed simulation would be skipped',
              pathlib.Path(paths[0] + '.completed').exists(), True)
        check('an unfinished one would not',
              pathlib.Path(paths[2] + '.completed').exists(), False)

        print('\nsetup_output_directory  (the real one, not a re-implementation)')
        exp_out = os.path.join(root, 'experiment_output')
        os.makedirs(exp_out, exist_ok=True)
        original_dir = om.dm.get_experiment_output_dir
        om.dm.get_experiment_output_dir = lambda name: exp_out
        seeded = os.path.join(root, 'cell', 'seed_0000.json')
        os.makedirs(os.path.dirname(seeded), exist_ok=True)
        try:
            os.environ.pop('SIMPLICITY_RESUME', None)
            out = om.setup_output_directory('exp', seeded)
            check('creates the directory on a first run',
                  os.path.isdir(out), True)
            pathlib.Path(os.path.join(out, 'stale.csv')).write_text('old')

            raised = ''
            try:
                om.setup_output_directory('exp', seeded)
            except RuntimeError as exc:
                raised = str(exc)
            check('refuses a second run by default',
                  raised.startswith('You already run an experiment'), True)
            check('and says how to resume', '--rerun' in raised, True)
            check('the earlier output is left alone',
                  sorted(os.listdir(out)), ['stale.csv'])

            os.environ['SIMPLICITY_RESUME'] = '1'
            out2 = om.setup_output_directory('exp', seeded)
            check('resuming returns the same directory', out2, out)
            check('and clears the abandoned attempt from it',
                  os.listdir(out2), [])
        finally:
            om.dm.get_experiment_output_dir = original_dir
            os.environ.pop('SIMPLICITY_RESUME', None)

        print('\nraise_on_stale_signals')
        pathlib.Path(paths[3] + '.started').touch()   # put a signal back
        blocked = False
        try:
            slurm.raise_on_stale_signals('irrelevant')
        except RuntimeError:
            blocked = True
        check('refuses to submit over an earlier attempt by default',
              blocked, True)
        os.environ['SIMPLICITY_RESUME'] = '1'
        blocked = False
        try:
            slurm.raise_on_stale_signals('irrelevant')
        except RuntimeError:
            blocked = True
        check('allows it when resuming', blocked, False)
    finally:
        slurm.sm.get_seeded_simulation_parameters_paths = original
        os.environ.pop('SIMPLICITY_RESUME', None)
        shutil.rmtree(root, ignore_errors=True)

    print('\n' + ('FAILURES:\n  ' + '\n  '.join(failures) if failures
                  else 'all passed'))
    sys.exit(1 if failures else 0)


if __name__ == '__main__':
    main()
