#!/usr/bin/env python3
'''reconcile_terminated_tasks must unblock a task Slurm reports as COMPLETED.

The case this covers, measured on profile_grid_#910: 120 tasks submitted, 120
profiles written and every one exited cleanly, `sacct` reporting COMPLETED for
all of them -- and only 116 completed signals on disk. The four without one
matched nothing in reconcile_terminated_tasks, because its only branch tested
SLURM_TERMINAL_FAILURE_STATES and COMPLETED is deliberately not in that set. So
poll_simulations_status' `left` never reached 0 and run_seeded_simulations
polled for hours, reporting a frozen progress snapshot against an empty queue.

Four things are asserted, the last two because getting them wrong would be
worse than the bug:

  a COMPLETED task with no state  -> completed, never failed. Slurm reports
                                     COMPLETED only after the job script exits
                                     0, so the simulation ran; marking it failed
                                     would under-report a successful run.
  a killed task                   -> failed, unchanged behaviour
  a task that recorded itself     -> untouched, not overwritten
  a task still RUNNING            -> untouched

    python tests/test_reconcile_completed.py
'''
import os
import shutil
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


def main():
    with Experiment(n_simulations=4) as exp:
        # four repeats, one per case, each in the STARTED state job() writes
        for position in range(4):
            exp.set(position, jobs.STARTED)
        exp.set(2, jobs.COMPLETED)          # recorded itself

        # the states sacct reports for each, in array-task order
        states = {'1': 'COMPLETED', '2': 'OUT_OF_MEMORY',
                  '3': 'COMPLETED', '4': 'RUNNING'}

        # Stub only what reaches outside this process: the job-id map job()
        # writes, and sacct itself. The repeats, their state and everything
        # between are the real thing on a real scratch Data/.
        original_index = slurm._build_slurm_id_map_index
        original_run = slurm.subprocess.run

        slurm._build_slurm_id_map_index = lambda name: {
            exp.key(position): ('4242', str(position + 1))
            for position in range(len(exp.pairs))}

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
            slurm.reconcile_terminated_tasks(exp.name, exp.pairs)
        finally:
            slurm._build_slurm_id_map_index = original_index
            slurm.subprocess.run = original_run

        print('\nreconcile_terminated_tasks against a real sacct reply')
        check('COMPLETED with no state is resolved as completed',
              exp.state(0), jobs.COMPLETED)
        check('OUT_OF_MEMORY is resolved as failed',
              exp.state(1), jobs.FAILED)
        check('a task that recorded itself is left alone',
              exp.state(2), jobs.COMPLETED)
        check('a RUNNING task is left alone', exp.state(3), jobs.STARTED)

        # a reconciled verdict says so, so a later reader can tell it apart
        # from a task that reported its own outcome
        group, index = exp.key(0)
        check('a reconciled repeat records why',
              jobs.get_state(exp.name, group, index).get('reconciled'),
              'slurm:COMPLETED')

        # The point of the fix: poll_simulations_status' `left` reaches 0.
        status = slurm.poll_simulations_status(exp.name, exp.pairs)
        check('repeats left unresolved (what the loop waits on)',
              status.left, 1)             # only the genuinely RUNNING one
        check('and the counts add up',
              (status.completed, status.failed), (2, 1))

    print('\n' + ('FAILURES:\n  ' + '\n  '.join(failures) if failures
                  else 'all passed'))
    sys.exit(1 if failures else 0)


if __name__ == '__main__':
    main()
