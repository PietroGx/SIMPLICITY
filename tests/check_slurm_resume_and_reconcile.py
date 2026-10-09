#!/usr/bin/env python3
'''Resume and reconciliation, against a real Slurm. Run this ON THE CLUSTER.

    python tests/check_slurm_resume_and_reconcile.py --exp-num 90
    python tests/check_slurm_resume_and_reconcile.py --exp-num 90 --skip-kill

WHY

Two paths in simplicity/runners/slurm.py are covered only by tests/fakeslurm:

  --rerun / resume    prepare_resume keeps completed repeats and clears the
                      rest; job() bails out before touching state for the ones
                      it keeps. The refactor changed what a resume is keyed on
                      -- a per-repeat state file instead of a .completed next
                      to a per-seed json -- and no real resume has run since.

  reconcile_terminated_tasks   a repeat the scheduler kills outright never
                      reaches job()'s except block, so it stays `started` and
                      the polling loop waits on it forever unless the
                      reconciler asks sacct and records the outcome. This is
                      the branch that unblocked profile_grid_#910. The clean
                      pipeline run never triggered it: nothing failed.

A fake can produce both on demand, which is why it does. What a fake cannot
confirm is that a real scheduler behaves the way the fake pretends -- and the
fake was already caught lying once, about squeue collapsing pending array
tasks into one row (fixed in v2.4.72).

WHAT IT DOES

Phase 1  submits four trivial repeats and lets them finish.
Phase 2  clears two of them and re-runs with SIMPLICITY_RESUME=1: the two that
         finished must be skipped, not recomputed, and the two cleared must run.
Phase 3  submits a second experiment whose repeats cannot finish inside a
         2-minute walltime, so Slurm kills them. The reconciler must mark them
         failed from sacct and the loop must exit instead of hanging.

         It is sized by WORK, not by horizon. The first version used a 3000-day
         final_time against 200 individuals and every task completed in twelve
         seconds -- final_time is simulated days, and a small population burns
         through them instantly. Nothing was killed, the reconciler never ran,
         and the phase passed anyway, because "no repeat left started" and
         "every repeat reached a terminal state" are both trivially true when
         everything completes. If nothing is killed this now reports
         INCONCLUSIVE and exits non-zero rather than going green.

Phase 3 deliberately produces FAILED repeats and burns a couple of minutes of
walltime on 4000 individuals. --skip-kill leaves it out.

It writes into Data/ under its own experiment names, so pick an --exp-num
nothing else is using. Everything it creates is left on disk for inspection.
'''
import argparse
import os
import shutil
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)

import simplicity.dir_manager as dm
import simplicity.jobs as jobs
import simplicity.settings_manager as sm

failures = []
inconclusive = []


def check(label, got, want):
    ok = got == want
    print(f"  [{'ok ' if ok else 'FAIL'}] {label}: {got!r}")
    if not ok:
        failures.append(f'{label}: got {got!r}, want {want!r}')
    return ok


def submit(name, final_time, seeds=2, walltime=None,
           population=200, infected=20):
    import simplicity.runme as runme
    import simplicity.runners.slurm as slurm
    if walltime:
        os.environ['SIMPLICITY_SLURM_TIME'] = walltime
    runme.run_experiment(
        name,
        lambda: ({'R': [1.1, 1.4]},
                 {'population_size': population,
                  'infected_individuals_at_start': infected,
                  'final_time': final_time}, seeds),
        simplicity_runner=slurm, archive_experiment=False)


def states(name):
    return {f'{g}/{r["index"]}': jobs.state_of(name, g, r['index'])
            for g, r in jobs.all_repeats(name)}


def output_mtimes(name):
    out = {}
    for sod in dm.get_simulation_output_dirs(name):
        for ssod in dm.get_seeded_simulation_output_dirs(sod):
            path = os.path.join(ssod, 'final_time.csv')
            if os.path.isfile(path):
                out[ssod] = os.path.getmtime(path)
    return out


def phase_1_and_2(exp_num):
    name = f'slurm_resume_check_#{exp_num}'
    print(f'\n=== phase 1: a clean run of four repeats  ({name})')
    submit(name, final_time=10)

    first = states(name)
    check('all four completed', sorted(set(first.values())), ['completed'])
    before = output_mtimes(name)
    check('all four wrote output', len(before), 4)

    print('\n=== phase 2: clear two, resume, and see what is recomputed')
    repeats = jobs.all_repeats(name)
    cleared = repeats[:2]
    kept = repeats[2:]
    for group, record in cleared:
        jobs.clear_state(name, group, record['index'])
        shutil.rmtree(os.path.join(dm.get_experiment_output_dir(name), group,
                                   record['stem'], jobs.repeat_label(record)),
                      ignore_errors=True)
    check('two repeats are now unstarted',
          sum(1 for v in states(name).values() if v is None), 2)

    time.sleep(1.1)        # so a recomputed file has a visibly newer mtime
    os.environ['SIMPLICITY_RESUME'] = '1'
    try:
        submit(name, final_time=10)
    finally:
        os.environ.pop('SIMPLICITY_RESUME', None)

    after = states(name)
    check('everything is completed again',
          sorted(set(after.values())), ['completed'])

    now = output_mtimes(name)
    recomputed = {p for p, t in now.items() if before.get(p, 0) < t - 0.5}
    expected = {os.path.join(dm.get_experiment_output_dir(name), g,
                             r['stem'], jobs.repeat_label(r))
                for g, r in cleared}
    check('exactly the cleared repeats were recomputed',
          recomputed == expected, True)
    if recomputed != expected:
        print(f'      recomputed: {sorted(os.path.basename(p) for p in recomputed)}')
        print(f'      expected  : {sorted(os.path.basename(p) for p in expected)}')
    check('and the kept ones were not touched',
          all(abs(now[p] - before[p]) < 0.5
              for g, r in kept
              for p in [os.path.join(dm.get_experiment_output_dir(name), g,
                                     r['stem'], jobs.repeat_label(r))]
              if p in now and p in before), True)


def phase_3(exp_num):
    name = f'slurm_kill_check_#{exp_num}'
    print(f'\n=== phase 3: tasks the scheduler kills  ({name})')
    print('    A 2-minute walltime against a population large enough that one')
    print('    repeat cannot finish inside it. Slurm then terminates the task,')
    print('    so job() never records a terminal state and only')
    print('    reconcile_terminated_tasks can unblock the loop.')
    print('    NOTE final_time is SIMULATED days, not wall time: a small')
    print('    population burns through a long horizon in seconds. The kill has')
    print('    to come from real work, which is why this one is big.')
    import simplicity.runners.slurm as slurm
    saved = slurm.RECONCILE_INTERVAL_S
    slurm.RECONCILE_INTERVAL_S = 60        # do not wait the default 15 minutes
    started = time.time()
    try:
        submit(name, final_time=1095, seeds=1, population=4000,
               infected=200, walltime='00:02:00')
    except Exception as exc:
        print(f'  [note] the run raised: {type(exc).__name__}: {exc}')
    finally:
        slurm.RECONCILE_INTERVAL_S = saved
        os.environ.pop('SIMPLICITY_SLURM_TIME', None)

    elapsed = time.time() - started
    final = states(name)
    print(f'    final states after {elapsed/60:.1f} min: {sorted(set(final.values()))}')
    check('the loop exited rather than hanging', elapsed < 60 * 25, True)
    check('no repeat was left started',
          [k for k, v in final.items() if v == jobs.STARTED], [])
    check('every repeat reached a terminal state',
          all(v in jobs.RESOLVED for v in final.values()), True)

    reconciled = [k for k, v in final.items() if v == jobs.FAILED
                  and (jobs.get_state(name, *_split(k)) or {}).get('reconciled')]
    print(f'    reconciled from sacct: {len(reconciled)} of {len(final)}')

    # If nothing was killed, the reconciler never ran and the three checks
    # above passed for the wrong reason -- they are all trivially true when
    # every task completes. Say so instead of going green: a test that reports
    # success while exercising nothing is worse than no test.
    if all(v == jobs.COMPLETED for v in final.values()):
        inconclusive.append(
            'phase 3 did not kill anything: every repeat completed inside the '
            'walltime, so reconcile_terminated_tasks was never reached. Raise '
            'the population or lower SIMPLICITY_SLURM_TIME and run it again.')
        print('    INCONCLUSIVE -- nothing was killed, so nothing was tested.')
        return

    check('a killed repeat was recorded failed',
          any(v == jobs.FAILED for v in final.values()), True)
    check('and the reconciler is what recorded it',
          len(reconciled) > 0, True)


def _split(key):
    group, _, index = key.rpartition('/')
    return group, int(index)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--exp-num', type=int, required=True)
    parser.add_argument('--skip-kill', action='store_true',
                        help='skip phase 3, which burns a minute of walltime '
                             'and leaves failed repeats behind')
    args = parser.parse_args()

    os.environ.setdefault('SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM', '4')
    os.environ.setdefault('SIMPLICITY_SLURM_MEM', '2G')
    os.environ.setdefault('SIMPLICITY_SLURM_TIME', '00:20:00')

    phase_1_and_2(args.exp_num)
    if not args.skip_kill:
        phase_3(args.exp_num)

    print('\n' + '=' * 70)
    if inconclusive:
        print('INCONCLUSIVE:')
        for line in inconclusive:
            print(f'  {line}')
    if failures:
        print(f'{len(failures)} FAILURE(S):')
        for line in failures:
            print(f'  {line}')
    elif not inconclusive:
        print('resume and reconciliation work against a real Slurm')
    sys.exit(1 if (failures or inconclusive) else 0)


if __name__ == '__main__':
    main()
