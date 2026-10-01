'''
Submitting over an earlier attempt's signal files must refuse, before sbatch.

On hpc-login02, cal_long #1 was submitted into a directory that still carried
.completed signals from an earlier attempt. poll_simulations_status reads those
files rather than Slurm, so the first poll reported
total=900 ... left=0 ... completed=900 before anything ran. The polling loop
never executed, release_simulations is only reached from inside it, and the
900-task array sat (JobHeldUser) forever while the pipeline walked on to

    RuntimeError: No valid long-calibration data found in ..._cal_long_#1

    python tests/test_stale_signals.py
'''
import os, pathlib, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import simplicity.runners.slurm as slurm

failures = []


def check(label, got, want):
    ok = got == want
    print(f"  [{'ok ' if ok else 'FAIL'}] {label}: {got!r}")
    if not ok:
        failures.append(f'{label}: got {got!r}, want {want!r}')


def seeded_paths(d, n=4):
    paths = []
    for i in range(n):
        p = os.path.join(d, f'seed_{i:04d}.json')
        pathlib.Path(p).touch()
        paths.append(p)
    return paths


def with_paths(paths):
    saved = slurm.sm.get_seeded_simulation_parameters_paths
    slurm.sm.get_seeded_simulation_parameters_paths = lambda exp: paths
    return saved


def restore(saved):
    slurm.sm.get_seeded_simulation_parameters_paths = saved


def test_clean_directory_passes():
    with tempfile.TemporaryDirectory() as d:
        saved = with_paths(seeded_paths(d))
        try:
            check('no signals found', slurm.find_stale_signals('probe'), {})
            slurm.raise_on_stale_signals('probe')   # must not raise
            check('a clean experiment submits', True, True)
        finally:
            restore(saved)


def test_completed_signals_refuse():
    '''The exact shape of the hpc-login02 failure.'''
    with tempfile.TemporaryDirectory() as d:
        paths = seeded_paths(d)
        for p in paths:
            for suffix in ('.submitted', '.released', '.started', '.completed'):
                pathlib.Path(p + suffix).touch()
        saved = with_paths(paths)
        try:
            counts = slurm.find_stale_signals('probe')
            check('every signal kind is counted', sorted(counts),
                  ['.completed', '.released', '.started', '.submitted'])
            check('counted once per seed', set(counts.values()), {len(paths)})
            raised = None
            try:
                slurm.raise_on_stale_signals('probe')
            except RuntimeError as e:
                raised = e
            check('it refuses', raised is not None, True)
            msg = str(raised)
            print('    ---\n    ' + msg.replace('\n', '\n    ') + '\n    ---')
            check('says nothing was submitted',
                  'Nothing has been submitted' in msg, True)
            check('names the counts', '4x .completed' in msg, True)
            check('offers a way out', '--exp-num' in msg and '-delete' in msg,
                  True)
        finally:
            restore(saved)


def test_a_single_leftover_signal_refuses():
    '''Even one .started miscounts the run, so the guard is not thresholded.'''
    with tempfile.TemporaryDirectory() as d:
        paths = seeded_paths(d)
        pathlib.Path(paths[2] + '.started').touch()
        saved = with_paths(paths)
        try:
            check('one leftover is found',
                  slurm.find_stale_signals('probe'), {'.started': 1})
            raised = None
            try:
                slurm.raise_on_stale_signals('probe')
            except RuntimeError as e:
                raised = e
            check('one leftover refuses', raised is not None, True)
        finally:
            restore(saved)


def test_guard_runs_before_sbatch():
    '''A refused run must leave no held array behind, so the check has to come
    before the subprocess call.'''
    with tempfile.TemporaryDirectory() as d:
        paths = seeded_paths(d)
        for p in paths:
            pathlib.Path(p + '.completed').touch()
        saved = with_paths(paths)
        saved_run = slurm.subprocess.run
        called = []
        try:
            slurm.subprocess.run = lambda *a, **k: called.append(a) or None
            raised = None
            try:
                slurm.submit_simulations('probe', test_guard_runs_before_sbatch,
                                         n=len(paths))
            except RuntimeError as e:
                raised = e
            check('submit_simulations refuses', raised is not None, True)
            check('sbatch was never called', called, [])
            leftover = [p for p in paths
                        if pathlib.Path(p + '.submitted').exists()]
            check('no new .submitted signals were written', leftover, [])
        finally:
            slurm.subprocess.run = saved_run
            restore(saved)


if __name__ == '__main__':
    for t in [test_clean_directory_passes, test_completed_signals_refuse,
              test_a_single_leftover_signal_refuses,
              test_guard_runs_before_sbatch]:
        print(f'\n{t.__name__}')
        t()
    print('\n' + ('FAILURES:\n  ' + '\n  '.join(failures) if failures
                  else 'all passed'))
    sys.exit(1 if failures else 0)
