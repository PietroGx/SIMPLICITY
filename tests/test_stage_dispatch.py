#!/usr/bin/env python3
'''Every pipeline stage reaches its dispatch, with the runner stubbed out.

    python tests/test_stage_dispatch.py

WHY

CLAUDE.md: "The config is correct" is not verification -- exercise the actual
function, stubbing the dispatch. v2.4.26 was a hotfix for a console printout
that still subscripted a renamed key and killed the stage before it ever
submitted anything.

It happened again. Stage 7a moved each sweep's axis inside its named group, so
the unbound cal_2 settings went from

    {'nucleotide_substitution_rate': [...]}
to
    {'_scenario_groups': [{'name': 'standard_only',
                           'nucleotide_substitution_rate': [...]}]}

and the progress banner four lines above the dispatch still read
varying['nucleotide_substitution_rate']. KeyError, before submission. Nobody
saw it because that stage had never been run.

A printout cannot be caught by reading a config or by checking that the groups
are right; it is only caught by CALLING main(). That is all this does: for
every stage, stub the runner and the analysis that needs real output, call
main(), and require that the dispatch was reached.

It also checks what was dispatched -- every group named, no name repeated, and
each group's swept axes non-empty -- because a stage that submits the wrong
grid fails just as expensively as one that does not submit at all, only later.

Nothing is simulated and nothing is submitted: the dispatch is replaced before
main() is called, so the recorded settings are the ones the real run would use.
'''
import argparse
import contextlib
import io
import os
import shutil
import sys
import tempfile

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'scripts'))
sys.path.insert(0, os.path.join(REPO, 'scripts', 'experiments'))

import simplicity.dir_manager as dm
import simplicity.settings_manager as sm

EXP_NUM = 9999
failures = []


def fail(label, detail):
    print(f'  [FAIL] {label}')
    for line in str(detail).splitlines():
        print(f'         {line}')
    failures.append(f'{label}: {detail}')


def ok(label, detail=''):
    print(f"  [ok ] {label}{': ' + detail if detail else ''}")


# --------------------------------------------------------------- the recorder
class Dispatch:
    '''Stands in for experiment_script_runner.run_experiment_phase_script and
    keeps what the stage asked to submit.

    The name matters. Until 8e the stages called run_experiment_script and
    this stubbed that; when they moved to run_experiment_phase_script the
    stub stopped matching and every stage RAN ITS SIMULATIONS FOR REAL --
    the suite still passed, taking 427s instead of 2s. A stub that silently
    stops intercepting is worse than no stub, so assert the attribute is
    there before replacing it.
    '''

    def __init__(self):
        self.calls = []

    def __call__(self, runner, exp_num, settings_func, experiment_name,
                 *extra, **kw):
        varying, fixed, n_seeds = settings_func()
        self.calls.append({'experiment': experiment_name,
                           'varying': varying, 'fixed': fixed,
                           'n_seeds': n_seeds})
        # what run_experiment_phase_script returns: this phase's groups
        return [g.get('name') for g in varying.get('_scenario_groups', [])]


def install(module, dispatch):
    '''Replace the module's dispatch, refusing if it does not have one.'''
    name = 'run_experiment_phase_script'
    if not hasattr(module, name):
        raise AssertionError(
            f'{module.__name__} has no {name} to stub -- the stage was '
            f'rewired and this test would have run real simulations.')
    setattr(module, name, dispatch)


def long_nsr_table():
    '''What cal_1 would have left behind: a calibrated NSR_long per cohort.'''
    import impact_long_shedders_config as config
    sp = sm.read_standard_parameters_values()
    return {pair: 5e-4
            for pair in config.unique_long_tau_r_long_pairs(sp)}


def frozen_rows(scenarios, sp, derive):
    '''What cal_2 would have left behind: one production row per scenario.'''
    rows = []
    for scenario in scenarios:
        frozen = derive(scenario, sp)
        rows.append({
            'scenario_name': frozen['scenario_name'],
            'tau_3_long': frozen['tau_3_long'],
            'long_shedders_ratio': frozen['long_shedders_ratio'],
            'susceptibility_long': frozen['susceptibility_long'],
            'R': 1.03, 'IH_virus_emergence_rate': 0.1,
            'R_long': frozen['R_long'] if frozen['R_long'] is not None else 1.1,
            'nucleotide_substitution_rate_long':
                5e-4 if frozen['long_shedders_ratio'] > 0 else float('nan'),
            'nucleotide_substitution_rate': 1.4e-4,
            'target_osr_long': 0.00205, 'model_type': 'exp',
            'long_calib_source': 'test', 'std_calib_source': 'test',
        })
    return rows


def write_frozen_table(setup_dir_template, filename, rows):
    directory = setup_dir_template.format(exp_num=EXP_NUM)
    os.makedirs(directory, exist_ok=True)
    pd.DataFrame(rows).to_csv(os.path.join(directory, filename), index=False)


# ------------------------------------------------------------------ the stages
def stage_cal_1(module_name, config_name):
    def run(dispatch):
        module = __import__(module_name)
        install(module, dispatch)
        module.plot_and_fit_long_nsr_calibration = lambda *a, **k: {}
        sys.argv = [module_name, '--exp-num', str(EXP_NUM),
                    '--runner', 'serial', '--seeds', '2']
        module.main()
    return run


def stage_cal_2():
    def run(dispatch):
        import impact_long_shedders_cal_2 as module
        install(module, dispatch)
        module.read_calibrated_long_nsr = lambda *a, **k: long_nsr_table()
        module.compute_standard_nsr = lambda *a, **k: 1.4e-4
        sys.argv = ['impact_long_shedders_cal_2',
                    '--exp-num', str(EXP_NUM), '--runner', 'serial',
                    '--seeds', '2', '--target-osr-long', '0.00205']
        module.main()
    return run


def stage_exp(module_name, config_name, extra_argv=()):      # noqa: ARG001
    def run(dispatch):
        config = __import__(config_name)
        module = __import__(module_name)
        sp = sm.read_standard_parameters_values()
        # take them from the module that actually reads the table
        write_frozen_table(module.SETUP_DIR_TEMPLATE, module.TABLE_FILENAME,
                           frozen_rows(config.SCENARIOS, sp,
                                       config.derive_scenario_params))
        install(module, dispatch)
        sys.argv = [module_name, '--exp-num', str(EXP_NUM),
                    '--runner', 'serial', '--seeds', '2', *extra_argv]
        module.main()
    return run


STAGES = [
    ('cal_1', stage_cal_1('impact_long_shedders_cal_1',
                          'impact_long_shedders_config'), 1),
    ('cal_2', stage_cal_2(), 1),
    # ONE dispatch since 8e, carrying all five scenarios as groups. It used
    # to be five, one experiment each, submitted from worker threads.
    ('exp  ', stage_exp('impact_long_shedders_exp',
                        'impact_long_shedders_config'), 1),
]


def inspect(label, call):
    '''What the stage asked to submit has to be submittable.'''
    varying = call['varying']
    if '_scenario_groups' not in varying:
        # a plain sweep: every axis is a varying parameter
        axes = {k: v for k, v in varying.items() if isinstance(v, list)}
        if not axes:
            fail(f'{label}: nothing varies', varying)
        return

    groups = varying['_scenario_groups']
    names = [g.get('name') for g in groups]
    if any(n is None for n in names):
        fail(f'{label}: a group has no name',
             'groups fall back to group_00.., and every reader that scopes by '
             f'name then has nothing to scope by: {names}')
    elif len(set(names)) != len(names):
        fail(f'{label}: two groups share a name', names)
    else:
        bad = [n for n in names if '/' in str(n) or os.sep in str(n)]
        if bad:
            fail(f'{label}: a group name would become a path', bad)
        else:
            ok(f'{label}: groups {names}')

    for group, name in zip(groups, names):
        axes = {k: v for k, v in group.items() if isinstance(v, list)}
        if any(len(v) == 0 for v in axes.values()):
            fail(f'{label}: group {name} has an empty sweep axis', axes)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--verbose', action='store_true',
                        help="show each stage's own console output")
    args = parser.parse_args()

    root = tempfile.mkdtemp(prefix='stage_dispatch_')
    data = os.path.join(root, 'Data')
    os.makedirs(data)
    saved_cwd = os.getcwd()
    # the stages build "Data/<...>_setup_data_#N" relative to the cwd
    os.chdir(root)
    dm.set_data_dir(data)
    sm._data_dir = data

    try:
        for label, run, expected_dispatches in STAGES:
            print(f'\n=== {label}')
            dispatch = Dispatch()
            buffer = io.StringIO()
            try:
                if args.verbose:
                    run(dispatch)
                else:
                    with contextlib.redirect_stdout(buffer):
                        run(dispatch)
            except SystemExit as exc:
                if exc.code:
                    fail(f'{label}: exited {exc.code} before finishing',
                         buffer.getvalue()[-600:])
                    continue
            except Exception as exc:                            # noqa: BLE001
                fail(f'{label}: {type(exc).__name__}: {exc}',
                     'last output before it died:\n' +
                     '\n'.join(buffer.getvalue().splitlines()[-6:]))
                continue

            if not dispatch.calls:
                fail(f'{label}: returned without dispatching anything',
                     buffer.getvalue()[-600:])
                continue

            if expected_dispatches is not None and \
                    len(dispatch.calls) != expected_dispatches:
                fail(f'{label}: dispatched {len(dispatch.calls)} experiment(s)',
                     f'expected {expected_dispatches}: '
                     f'{[c["experiment"] for c in dispatch.calls]}')
                continue

            ok(f'{label}: dispatched',
               ', '.join(c['experiment'] for c in dispatch.calls))
            for call in dispatch.calls:
                inspect(f"{label} {call['experiment']}", call)
    finally:
        os.chdir(saved_cwd)
        shutil.rmtree(root, ignore_errors=True)

    print('\n' + '=' * 70)
    if failures:
        print(f'{len(failures)} FAILURE(S):')
        for line in failures:
            print(f'  {line}')
    else:
        print('every stage reaches its dispatch and submits named, '
              'non-empty groups')
    sys.exit(1 if failures else 0)


if __name__ == '__main__':
    main()
