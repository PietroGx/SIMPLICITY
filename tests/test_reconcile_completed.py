#!/usr/bin/env python3
'''reconcile_terminated_tasks must unblock a task Slurm reports as COMPLETED.

The case this covers, measured on profile_grid_#910: 120 tasks submitted, 120
profiles written and every one exited cleanly, `sacct` reporting COMPLETED for
all of them -- and only 116 .completed signals on disk. The four without one
matched nothing in reconcile_terminated_tasks, because its only branch tested
SLURM_TERMINAL_FAILURE_STATES and COMPLETED is deliberately not in that set. So
poll_simulations_status' `left` never reached 0 and run_seeded_simulations
polled for hours, reporting a frozen .progress snapshot against an empty queue.

Four things are asserted, the last two because getting them wrong would be
worse than the bug:

  a COMPLETED task with no signal      -> .completed, never .failed. Slurm
                                          reports COMPLETED only after the job
                                          script exits 0, so the simulation ran;
                                          marking it .failed would under-report
                                          a successful run.
  a killed task                        -> .failed, unchanged behaviour
  a task that signaled itself          -> untouched, no second signal
  a task still RUNNING                 -> untouched

    python tests/test_reconcile_completed.py
'''
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import simplicity.runners.slurm as slurm

failures = []


def check(label, got, want):
    ok = got == want
    print(f"  [{'ok ' if ok else 'FAIL'}] {label}: {got!r}")
    if not ok:
        failures.append(f'{label}: got {got!r}, want {want!r}')


def signals(path):
    return {suffix for suffix in ('.started', '.completed', '.failed')
            if pathlib.Path(path + suffix).exists()}


def main():
    root = tempfile.mkdtemp(prefix='reconcile_test_')
    try:
        # four tasks, one per case, each with the .started signal job() writes
        names = ['completed_no_signal', 'killed', 'signaled_itself', 'running']
        paths = [os.path.join(root, f'{name}.json') for name in names]
        for path in paths:
            pathlib.Path(path).write_text('{}')
            pathlib.Path(path + '.started').touch()
        pathlib.Path(paths[2] + '.completed').touch()   # signaled itself

        # the states sacct reports for each, in array-task order
        states = {'1': 'COMPLETED', '2': 'OUT_OF_MEMORY',
                  '3': 'COMPLETED', '4': 'RUNNING'}

        # Stub the three things that reach outside this process: the parameter
        # paths, the job-id map job() writes, and sacct itself. Everything
        # between them is the real function.
        original_paths = slurm.sm.get_seeded_simulation_parameters_paths
        original_index = slurm._build_slurm_id_map_index
        original_run = slurm.subprocess.run

        slurm.sm.get_seeded_simulation_parameters_paths = lambda name: paths
        slurm._build_slurm_id_map_index = lambda name: {
            path: ('4242', str(i + 1)) for i, path in enumerate(paths)}

        class Result:
            returncode = 0
            stdout = '\n'.join(f'4242_{task}|{state}'
                               for task, state in states.items()).encode()

        def fake_run(args, *rest, **kwargs):
            if args and 'sacct' in str(args[0]):
                return Result()
            return original_run(args, *rest, **kwargs)

        slurm.subprocess.run = fake_run
        try:
            slurm.reconcile_terminated_tasks('irrelevant')
        finally:
            slurm.sm.get_seeded_simulation_parameters_paths = original_paths
            slurm._build_slurm_id_map_index = original_index
            slurm.subprocess.run = original_run

        print('\nreconcile_terminated_tasks against a real sacct reply')
        check('COMPLETED with no signal is resolved as completed',
              signals(paths[0]), {'.started', '.completed'})
        check('OUT_OF_MEMORY is resolved as failed',
              signals(paths[1]), {'.started', '.failed'})
        check('a task that signaled itself is left alone',
              signals(paths[2]), {'.started', '.completed'})
        check('a RUNNING task is left alone',
              signals(paths[3]), {'.started'})

        # The point of the fix: poll_simulations_status' `left` reaches 0.
        resolved = sum(1 for path in paths
                       if signals(path) & {'.completed', '.failed'})
        check('tasks left unresolved (what the loop waits on)',
              len(paths) - resolved, 1)   # only the genuinely RUNNING one
    finally:
        shutil.rmtree(root, ignore_errors=True)

    print('\n' + ('FAILURES:\n  ' + '\n  '.join(failures) if failures
                  else 'all passed'))
    sys.exit(1 if failures else 0)


if __name__ == '__main__':
    main()
