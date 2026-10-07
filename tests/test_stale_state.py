#!/usr/bin/env python3
'''Submitting over an earlier attempt's repeat state must refuse, before sbatch.

On hpc-login02, cal_long #1 was submitted into a directory that still carried
completed signals from an earlier attempt. poll_simulations_status reads those
files rather than Slurm, so the first poll reported
total=900 ... left=0 ... completed=900 before anything ran. The polling loop
never executed, release_simulations is only reached from inside it, and the
900-task array sat (JobHeldUser) forever while the pipeline walked on to

    RuntimeError: No valid long-calibration data found in ..._cal_long_#1

Was tests/test_stale_signals.py, against five signal files per simulation. The
state now lives in one file per repeat, so the guard counts states instead of
suffixes -- but the failure it prevents is the same one.

    python tests/test_stale_state.py
'''
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

import simplicity.jobs as jobs
import simplicity.runners.slurm as slurm
from _experiment_fixture import Experiment

failures = []


def check(label, got, want):
    ok = got == want
    print(f"  [{'ok ' if ok else 'FAIL'}] {label}: {got!r}")
    if not ok:
        failures.append(f'{label}: got {got!r}, want {want!r}')


def test_clean_experiment_passes():
    with Experiment() as exp:
        check('no state found', slurm.find_stale_state(exp.name), {})
        slurm.raise_on_stale_state(exp.name)       # must not raise
        check('a clean experiment submits', True, True)


def test_completed_state_refuses():
    '''The exact shape of the hpc-login02 failure.'''
    with Experiment() as exp:
        for position in range(len(exp.pairs)):
            exp.set(position, jobs.COMPLETED)

        counts = slurm.find_stale_state(exp.name)
        check('the completed repeats are counted',
              counts, {jobs.COMPLETED: len(exp.pairs)})

        raised = None
        try:
            slurm.raise_on_stale_state(exp.name)
        except RuntimeError as exc:
            raised = exc
        check('it refuses', raised is not None, True)
        message = str(raised)
        print('    ---\n    ' + message.replace('\n', '\n    ') + '\n    ---')
        check('says nothing was submitted',
              'Nothing has been submitted' in message, True)
        check('names the counts', f'{len(exp.pairs)}x completed' in message, True)
        check('offers a way out',
              '--exp-num' in message and 'rm -rf' in message, True)


def test_a_single_leftover_refuses():
    '''Even one started repeat miscounts the run, so the guard is not
    thresholded.'''
    with Experiment() as exp:
        exp.set(2, jobs.STARTED)
        check('one leftover is found',
              slurm.find_stale_state(exp.name), {jobs.STARTED: 1})
        raised = None
        try:
            slurm.raise_on_stale_state(exp.name)
        except RuntimeError as exc:
            raised = exc
        check('one leftover refuses', raised is not None, True)


def test_guard_runs_before_sbatch():
    '''A refused run must leave no held array AND no submission map behind, so
    the check has to come before both.'''
    with Experiment() as exp:
        for position in range(len(exp.pairs)):
            exp.set(position, jobs.COMPLETED)

        original_run = slurm.subprocess.run
        called = []
        try:
            slurm.subprocess.run = lambda *a, **k: called.append(a) or None
            raised = None
            try:
                slurm.run_seeded_simulations(exp.name,
                                             test_guard_runs_before_sbatch)
            except RuntimeError as exc:
                raised = exc
            check('run_seeded_simulations refuses', raised is not None, True)
            check('sbatch was never called', called, [])
            check('no submission map was written',
                  os.listdir(jobs.submissions_dir(exp.name)), [])
            check('no repeat was moved to submitted',
                  [s for s in exp.states() if s == jobs.SUBMITTED], [])
        finally:
            slurm.subprocess.run = original_run


def test_resume_does_not_refuse():
    '''A resume is the one case where earlier state is the point.'''
    with Experiment() as exp:
        exp.set(0, jobs.COMPLETED)
        exp.set(1, jobs.FAILED)
        exp.set(2, jobs.STARTED)
        os.environ['SIMPLICITY_RESUME'] = '1'
        try:
            slurm.raise_on_stale_state(exp.name)   # must not raise
            check('resuming is allowed over earlier state', True, True)
            check('the completed one is kept', exp.state(0), jobs.COMPLETED)
            check('the failed one is cleared', exp.state(1), None)
            check('the half-started one is cleared', exp.state(2), None)
        finally:
            os.environ.pop('SIMPLICITY_RESUME', None)


if __name__ == '__main__':
    for test in [test_clean_experiment_passes, test_completed_state_refuses,
                 test_a_single_leftover_refuses, test_guard_runs_before_sbatch,
                 test_resume_does_not_refuse]:
        print(f'\n{test.__name__}')
        test()
    print('\n' + ('FAILURES:\n  ' + '\n  '.join(failures) if failures
                  else 'all passed'))
    sys.exit(1 if failures else 0)
