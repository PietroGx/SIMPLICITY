#!/usr/bin/env python3
'''Resuming an experiment must keep what finished and re-run only the rest.

Without a resume, re-submitting an experiment is refused twice over:
raise_on_stale_state stops it before sbatch, and setup_output_directory raises
"You already run an experiment with the same name!" per repeat. Both are right
by default -- which of two runs matters is not something the code can know.
SIMPLICITY_RESUME=1 (the --rerun flag) says the earlier one matters and the gap
should be filled.

What that has to get right, and what this checks:

  completed kept      a repeat recorded COMPLETED is left alone and its output
                      is never recomputed. Cell #17 of the calibrated grid
                      stalled with 61 of 150 finished; throwing those away was
                      the first thing I built and it was wrong.

  failed retried      a FAILED repeat is work still to do on a resume, not a
                      settled answer, so its state is cleared too.

  partial cleared     a repeat that started and never finished has a
                      half-written output directory. It is removed before the
                      re-run: leaving it mixes the abandoned attempt's files
                      with the new ones and nothing downstream can tell which
                      is which.

  progress cleared    the stale progress snapshot goes with the state, or the
                      monitor reports hours of elapsed time for a repeat that
                      has only just restarted.

  refuses by default  with SIMPLICITY_RESUME unset, both guards still fire.
                      A resume must be asked for, never inferred.

    python tests/test_resume.py
'''
import os
import pathlib
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

import simplicity.jobs as jobs
import simplicity.output_manager as om
import simplicity.runners.slurm as slurm
from _experiment_fixture import Experiment

failures = []


def check(label, got, want):
    ok = got == want
    print(f"  [{'ok ' if ok else 'FAIL'}] {label}: {got!r}")
    if not ok:
        failures.append(f'{label}: got {got!r}, want {want!r}')


def main():
    test_renumbering_under_state_is_refused()
    try:
        with Experiment(n_simulations=4) as exp:
            # one of each: finished, failed, half-started, never touched
            exp.set(0, jobs.COMPLETED)
            exp.set(1, jobs.FAILED)
            exp.set(2, jobs.STARTED)
            # the half-started one also left a progress snapshot behind
            group2, index2 = exp.key(2)
            pathlib.Path(jobs.progress_path(exp.name, group2, index2)
                         ).write_text('{"time": 400}')

            print('\nprepare_resume')
            kept, to_run = slurm.prepare_resume(exp.name)
            check('one finished repeat is kept', kept, 1)
            check('the other three are to run', to_run, 3)
            check('the completed one is untouched',
                  exp.state(0), jobs.COMPLETED)
            check('the failed one is cleared back to unstarted',
                  exp.state(1), None)
            check('the half-started one is cleared', exp.state(2), None)
            check('the untouched one is still untouched', exp.state(3), None)
            check('and its stale progress snapshot is gone',
                  jobs.get_progress(exp.name, group2, index2), None)

            print('\nslurm.job skips what finished')
            os.environ['SIMPLICITY_RESUME'] = '1'
            check('resuming() reads the environment', slurm.resuming(), True)
            group0, index0 = exp.key(0)
            check('a completed repeat would be skipped',
                  jobs.state_of(exp.name, group0, index0), jobs.COMPLETED)
            check('an unfinished one would not',
                  jobs.state_of(exp.name, group2, index2), None)
            os.environ.pop('SIMPLICITY_RESUME', None)

            print('\nsetup_output_directory  (the real one)')
            group, record = jobs.all_repeats(exp.name)[0]
            out = om.setup_output_directory(exp.name, group, record)
            check('creates the directory on a first run',
                  os.path.isdir(out), True)
            check('under the group, named for the simulation and the seed',
                  (os.path.basename(out),
                   os.path.basename(os.path.dirname(os.path.dirname(out)))),
                  (jobs.repeat_label(record), group))
            pathlib.Path(os.path.join(out, 'stale.csv')).write_text('old')

            raised = ''
            try:
                om.setup_output_directory(exp.name, group, record)
            except RuntimeError as exc:
                raised = str(exc)
            check('refuses a second run by default',
                  raised.startswith('You already run an experiment'), True)
            check('and says how to resume', '--rerun' in raised, True)
            check('the earlier output is left alone',
                  sorted(os.listdir(out)), ['stale.csv'])

            os.environ['SIMPLICITY_RESUME'] = '1'
            out2 = om.setup_output_directory(exp.name, group, record)
            check('resuming returns the same directory', out2, out)
            check('and clears the abandoned attempt from it',
                  os.listdir(out2), [])
            # the simulation's own parameters.json sits one level up, written at
            # setup -- the rmtree must not take it with the seed directory
            check("the simulation's parameters.json survives",
                  os.path.isfile(os.path.join(os.path.dirname(out),
                                              'parameters.json')), True)
            os.environ.pop('SIMPLICITY_RESUME', None)

            print('\nraise_on_stale_state')
            exp.set(3, jobs.STARTED)        # put some state back
            blocked = False
            try:
                slurm.raise_on_stale_state(exp.name)
            except RuntimeError:
                blocked = True
            check('refuses to submit over an earlier attempt by default',
                  blocked, True)
            os.environ['SIMPLICITY_RESUME'] = '1'
            blocked = False
            try:
                slurm.raise_on_stale_state(exp.name)
            except RuntimeError:
                blocked = True
            check('allows it when resuming', blocked, False)
    finally:
        os.environ.pop('SIMPLICITY_RESUME', None)

    print('\n' + ('FAILURES:\n  ' + '\n  '.join(failures) if failures
                  else 'all passed'))
    sys.exit(1 if failures else 0)


def test_renumbering_under_state_is_refused():
    '''Changing n_seeds renumbers repeats, so resuming would mis-attribute.'''
    import simplicity.settings_manager as sm
    print('\nwrite_repeats refuses to renumber a group that carries state')
    with Experiment(n_simulations=3, n_seeds=2) as exp:
        before = jobs.repeat_identity(jobs.read_repeats(exp.name, 'main'))
        exp.set(0, jobs.COMPLETED)

        # the same n_seeds is not a renumbering, so a resume still works
        jobs.write_repeats(exp.name)
        check('rewriting the same list is allowed',
              jobs.repeat_identity(jobs.read_repeats(exp.name, 'main')), before)
        check('and the completed repeat is untouched',
              exp.state(0), jobs.COMPLETED)

        # now change n_seeds, which moves everything after simulation 0
        record = sm.read_experiment_settings_file(exp.name)
        record['n_seeds'] = 4
        for group in record['groups']:
            group['n_seeds'] = 4
        with open(sm.get_experiment_settings_file_path(exp.name), 'w') as handle:
            import json; json.dump(record, handle)

        raised = ''
        try:
            jobs.write_repeats(exp.name)
        except RuntimeError as exc:
            raised = str(exc)
        check('renumbering under live state is refused', bool(raised), True)
        check('and it says why', 'WRONG repeats' in raised, True)
        check('and how to proceed', '--exp-num' in raised, True)
        check('the repeat list on disk is unchanged',
              jobs.repeat_identity(jobs.read_repeats(exp.name, 'main')), before)

        # with the state cleared, the same change is fine
        jobs.clear_state(exp.name, 'main', 0)
        jobs.write_repeats(exp.name)
        check('clearing the state allows it',
              len(jobs.read_repeats(exp.name, 'main')), 12)


if __name__ == '__main__':
    main()
