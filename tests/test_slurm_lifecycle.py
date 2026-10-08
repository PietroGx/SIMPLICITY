#!/usr/bin/env python3
'''The whole Slurm control plane, end to end, against a fake controller.

Not a unit test. simplicity.runners.slurm really submits, polls, releases,
reconciles and terminates; the tasks really run as subprocesses and really
resolve their own repeat. Only sbatch/squeue/scontrol/sacct are fake -- see
tests/fakeslurm/_fakeslurm.py for why a fake beats a real local controller for
the cases that matter.

WHAT THIS IS FOR

The refactor changed how a task finds its work. It used to index an unsorted
os.walk, performed separately by the submitting process and by each task. It now
looks its array position up in a map the submitter recorded. That is strictly
better, and it is also the single place where a mistake would quietly run the
wrong parameters under the right folder name -- the failure this repo has burned
four production runs on. Nothing short of really submitting an array and really
running every task checks it.

So the first and most important assertion here is dull: every array position
resolved to a distinct repeat, every repeat ran exactly once, and none ran twice.

Then the four externally-killed cases, each of which has hung a real run:

    hang         SIGKILLed mid-work -> stays `started`, sacct says TIMEOUT,
                 reconcile_terminated_tasks must mark it failed
    silent       exits 0 but its terminal state write is lost, sacct says
                 COMPLETED -> must be marked completed, never failed
                 (profile_grid_#910: 116 of 120)
    raise        the work raises -> job() records failed itself
    launch_fail  Slurm requeues it held -> re-released, then runs

    python tests/test_slurm_lifecycle.py
'''
import json
import os
import shutil
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, REPO)

import simplicity.dir_manager as dm
import simplicity.jobs as jobs
import simplicity.settings_manager as sm
import simplicity.runners.slurm as slurm

FAKE_BIN = os.path.join(HERE, 'fakeslurm')
WORK = '_fake_work.run_seeded_simulation'

failures = []


def check(label, got, want):
    ok = got == want
    print(f"  [{'ok ' if ok else 'FAIL'}] {label}: {got!r}")
    if not ok:
        failures.append(f'{label}: got {got!r}, want {want!r}')


def build(data_dir, name, n_simulations, n_seeds):
    dm.set_data_dir(data_dir)
    sm._data_dir = data_dir
    dm.create_directories(name)
    varying = {'nucleotide_substitution_rate':
               [1e-4 * (i + 1) for i in range(n_simulations)]}
    sm.write_experiment_settings(
        name,
        sm.generate_experiment_settings(varying, {'population_size': 20,
                                                  'final_time': 2}),
        n_seeds)
    sm.read_settings_and_write_simulation_parameters(name)
    jobs.write_repeats(name)


def drive(root, name, plan, cap=2, kill_state='TIMEOUT'):
    """Run the real polling loop against the fake controller."""
    data_dir = os.path.join(root, 'Data')
    ran_dir = os.path.join(root, 'ran')
    state_file = os.path.join(root, 'fakeslurm.json')

    saved_env = dict(os.environ)
    saved_cwd = os.getcwd()
    saved = (slurm.RECONCILE_INTERVAL_S,
             slurm.LAUNCH_FAILURE_RECONCILE_INTERVAL_S,
             slurm.LONG_RUNNING_REPORT_INTERVAL_S,
             time.sleep)
    try:
        # cwd decides the data directory: dir_manager computes it from os.getcwd
        # at import, so the task subprocesses find the same scratch Data/ without
        # being told about it.
        os.chdir(root)
        os.environ['PATH'] = FAKE_BIN + os.pathsep + os.environ['PATH']
        os.environ['PYTHONPATH'] = os.pathsep.join([REPO, HERE])
        os.environ['FAKESLURM_STATE'] = state_file
        os.environ['FAKESLURM_PLAN'] = json.dumps(plan)
        os.environ['FAKESLURM_KILL_STATE'] = kill_state
        os.environ['FAKESLURM_REPO'] = REPO
        os.environ['FAKESLURM_DATA_DIR'] = data_dir
        os.environ['FAKEWORK_RAN_DIR'] = ran_dir
        os.environ['SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM'] = str(cap)
        os.environ['SIMPLICITY_SLURM_MEM'] = '1G'
        os.environ['SIMPLICITY_SLURM_TIME'] = '00:10:00'
        os.environ.pop('SIMPLICITY_RESUME', None)

        # reconcile every pass, and do not wait seven seconds between them: both
        # are tuning constants, not logic
        slurm.RECONCILE_INTERVAL_S = 0
        slurm.LAUNCH_FAILURE_RECONCILE_INTERVAL_S = 0
        slurm.LONG_RUNNING_REPORT_INTERVAL_S = 10 ** 9
        time.sleep = lambda _seconds: saved[3](0.01)

        import _fake_work
        slurm.run_seeded_simulations(name, _fake_work.run_seeded_simulation)
    finally:
        (slurm.RECONCILE_INTERVAL_S,
         slurm.LAUNCH_FAILURE_RECONCILE_INTERVAL_S,
         slurm.LONG_RUNNING_REPORT_INTERVAL_S,
         time.sleep) = saved
        os.chdir(saved_cwd)
        os.environ.clear()
        os.environ.update(saved_env)


def ran(root):
    """{group_index: times run}, from the markers the fake work writes."""
    path = os.path.join(root, 'ran')
    out = {}
    for name in sorted(os.listdir(path)) if os.path.isdir(path) else []:
        with open(os.path.join(path, name)) as handle:
            out[name] = len([line for line in handle if line.strip()])
    return out


# ---------------------------------------------------------------------- tests

def test_every_position_resolves_to_its_own_repeat():
    '''The assertion the refactor lives or dies on.'''
    root = tempfile.mkdtemp(prefix='slurm_e2e_')
    try:
        build(os.path.join(root, 'Data'), 'probe', n_simulations=3, n_seeds=2)
        expected = [(group, record['index'])
                    for group, record in jobs.all_repeats('probe')]
        drive(root, 'probe', plan={})

        print('\nevery array position resolved to its own repeat')
        check('six repeats were defined', len(expected), 6)
        markers = ran(root)
        check('every repeat ran',
              sorted(markers), sorted(f'{g}_{i:06d}' for g, i in expected))
        check('and none ran more than once',
              sorted(set(markers.values())), [1])

        print('\nand every one is recorded completed')
        check('states', jobs.count_states('probe'), {jobs.COMPLETED: 6})

        print('\nthe repeat order is the settings order, not readdir order')
        check('array position n ran simulation n//n_seeds, seed n%n_seeds',
              [(r['simulation'], r['seed'])
               for _g, r in jobs.all_repeats('probe')],
              [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)])
    finally:
        shutil.rmtree(root, ignore_errors=True)


def test_externally_killed_and_lost_writes_are_reconciled():
    root = tempfile.mkdtemp(prefix='slurm_e2e_')
    try:
        build(os.path.join(root, 'Data'), 'probe', n_simulations=4, n_seeds=1)
        # array positions are 1-based: task 1..4 -> repeat index 0..3
        plan = {'1': 'ok', '2': 'hang', '3': 'silent', '4': 'raise'}
        drive(root, 'probe', plan=plan, kill_state='OUT_OF_MEMORY')

        print('\neach planned fate reached its terminal state')
        states = [jobs.state_of('probe', 'main', i) for i in range(4)]
        check('ok          -> completed', states[0], jobs.COMPLETED)
        check('hang        -> failed, reconciled from sacct', states[1],
              jobs.FAILED)
        check('silent      -> completed, reconciled from sacct', states[2],
              jobs.COMPLETED)
        check('raise       -> failed, recorded by job() itself', states[3],
              jobs.FAILED)

        print('\nthe two reconciled ones say so; the two self-reported do not')
        reconciled = [bool((jobs.get_state('probe', 'main', i) or {})
                           .get('reconciled')) for i in range(4)]
        check('reconciled flags', reconciled, [False, True, True, False])

        print('\nthe loop terminated, which is the whole point')
        status = slurm.poll_simulations_status(
            'probe', [['main', i] for i in range(4)])
        check('nothing left unresolved', status.left, 0)
        check('counts', (status.completed, status.failed), (2, 2))

        # A killed task must never be recorded completed, and an exit-0 task
        # must never be recorded failed: the first loses a run, the second
        # fabricates one.
        print('\nthe two directions are not interchangeable')
        check('the killed task was not called completed',
              states[1] == jobs.COMPLETED, False)
        check('the exit-0 task was not called failed',
              states[2] == jobs.FAILED, False)
    finally:
        shutil.rmtree(root, ignore_errors=True)


def test_launch_failures_are_re_released():
    """TWO tasks fail to launch, so they share a reason and squeue collapses
    them into one row with a range in the task field.

    One failing task can report a bare integer, which the old int() parse
    happened to accept. Two make the row read "2-3", and every pending task
    was then silently skipped -- the repeat sat held and the polling loop
    waited on it forever. Found by check_slurm_interface.py on slurm 26.05.4,
    not here: these shims used to emit one row per task.
    """
    root = tempfile.mkdtemp(prefix='slurm_e2e_')
    try:
        build(os.path.join(root, 'Data'), 'probe', n_simulations=3, n_seeds=1)
        drive(root, 'probe', plan={'2': 'launch_fail', '3': 'launch_fail'})

        print('\ntasks Slurm requeued held are re-released, then run')
        check('all three repeats ended completed',
              [jobs.state_of('probe', 'main', i) for i in range(3)],
              [jobs.COMPLETED] * 3)
        for index in (1, 2):
            record = jobs.get_state('probe', 'main', index) or {}
            check(f'repeat {index}: the re-releases were counted',
                  record.get('attempts', 0) >= 1, True)
        check('and they really did run', sorted(ran(root)),
              ['main_000000', 'main_000001', 'main_000002'])
    finally:
        shutil.rmtree(root, ignore_errors=True)


if __name__ == '__main__':
    for test in [test_every_position_resolves_to_its_own_repeat,
                 test_externally_killed_and_lost_writes_are_reconciled,
                 test_launch_failures_are_re_released]:
        print(f'\n{"=" * 70}\n{test.__name__}')
        test()
    print('\n' + ('FAILURES:\n  ' + '\n  '.join(failures) if failures
                  else 'all passed'))
    sys.exit(1 if failures else 0)
