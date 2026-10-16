#!/usr/bin/env python3
'''The simulation's graceful stop must sit just under the Slurm wall.

    python tests/test_runtime_guard.py

WHY

extrande's max_runtime (extrande.py:259) ENDS THE RUN AND WRITES ITS OUTPUT.
A Slurm walltime kill does not: the task dies mid-write and the repeat has to
be reconciled and re-run. So whichever of the two is lower decides what
happens to a long simulation, and the useful ordering is max_runtime first.

It was 86000s (~24h), chosen against a 1-day wall. The wall then went to 2
days (the guard binding again at ~24h), and in v2.4.100/101 to 5h -- at which
point the guard sat FOUR AND A HALF HOURS ABOVE the wall and could never fire.
Nothing announced that; two numbers set in two files drifted past each other.

So the guard is no longer a number. It is derived from SIMPLICITY_SLURM_TIME
by dir_manager.graceful_max_runtime(), and every builder STATES it in
fixed_params rather than inheriting it -- because read_standard_parameters_
values lets a machine-local standard_values.json win over the code default,
so a code-only change would not reach an existing install. That is the same
trap v2.4.93 hit with the consensus mode.

Checked here: the parse, the margin, and that all three stages really write a
value that follows --slurm-time.
'''
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'scripts'))
sys.path.insert(0, os.path.join(REPO, 'scripts', 'experiments'))

import pandas as pd

import simplicity.dir_manager as dm
import simplicity.settings_manager as sm
import impact_long_shedders_config as cfg

failures = []


def check(label, got, want):
    ok = got == want
    print(f"  [{'ok ' if ok else 'FAIL'}] {label}: {got!r}")
    if not ok:
        failures.append(f'{label}: got {got!r}, want {want!r}')


# ===========================================================================
def test_parse():
    print('\n=== every walltime form sbatch accepts')
    for spec, want in (('05:00:00', 18000), ('12:00:00', 43200),
                       ('2-00:00:00', 172800), ('1-12:00:00', 129600),
                       ('00:20:00', 1200), ('30:00', 1800), ('45', 2700),
                       ('3-06', 280800), ('0-00:00:30', 30)):
        check(f'{spec} -> seconds', dm.slurm_seconds(spec), want)


def test_margin():
    print('\n=== the guard lands under the wall, with room to write')
    for spec in ('05:00:00', '12:00:00', '2-00:00:00', '00:20:00'):
        wall = dm.slurm_seconds(spec)
        guard = dm.graceful_max_runtime(spec)
        check(f'{spec}: guard is BELOW the wall', guard < wall, True)
        check(f'{spec}: margin >= 10 min or 5%',
              wall - guard >= max(dm.RUNTIME_MARGIN_S,
                                  int(wall * dm.RUNTIME_MARGIN_FRACTION)), True)
    # the live default, stated so a change to it is visible here
    check('at the 5h default the guard is 4h45m',
          dm.graceful_max_runtime('05:00:00'), 17100)
    check('and a tiny wall still leaves a usable guard',
          dm.graceful_max_runtime('00:01:00') >= 60, True)


# ===========================================================================
def frozen_rows():
    sp = sm.read_standard_parameters_values()
    rows = []
    for scenario in cfg.SCENARIOS:
        frozen = cfg.derive_scenario_params(scenario, sp)
        rows.append(pd.Series({
            **frozen, 'R': 1.03, 'IH_virus_emergence_rate': 0.1,
            'nucleotide_substitution_rate': 1.4e-4,
            'nucleotide_substitution_rate_long':
                5e-4 if frozen['is_long'] else float('nan')}))
    return rows


RANGES = {'min': 1e-4, 'max': 1e-2, 'steps': 3}


def builders():
    """Every stage that submits simulations, by the name its script uses."""
    return {
        'cal_1': cfg.build_cal1_settings(2, RANGES),
        'cal_2': cfg.build_cal2_settings(2, RANGES),
        'exp':   cfg.build_exp_phase_settings(frozen_rows(), 3),
    }


def test_every_stage_states_it():
    print('\n=== all three stages write the guard into fixed_params')
    for stage, build in builders().items():
        _, fixed, _ = build()
        check(f'{stage} states max_runtime',
              'max_runtime' in fixed, True)
        check(f'{stage} is under the live wall',
              fixed.get('max_runtime') < dm.slurm_seconds(
                  os.environ['SIMPLICITY_SLURM_TIME']), True)


def test_it_reaches_a_simulation():
    """fixed_params is not the deliverable -- a simulation's own parameter
    record is. Every stage goes through generate_experiment_settings here,
    which is what the output directories are built from."""
    print('\n=== and every simulation record carries it')
    wall = dm.slurm_seconds(os.environ['SIMPLICITY_SLURM_TIME'])
    for stage, build in builders().items():
        varying, fixed, _ = build()
        records = sm.generate_experiment_settings(varying, fixed)
        values = {r['parameters'].get('max_runtime') for r in records}
        check(f'{stage}: every one of its {len(records)} simulations has one',
              values == {dm.graceful_max_runtime()}, True)
        check(f'{stage}: and it is under the wall',
              all(v is not None and v < wall for v in values), True)


def test_it_follows_the_flag():
    '''--slurm-time reaches the env through set_slurm_resource_env, and must
    reach the guard WITH it -- including for a builder constructed first,
    which is the order impact_long_shedders_exp uses.'''
    print('\n=== the guard follows --slurm-time')
    before = os.environ['SIMPLICITY_SLURM_TIME']
    try:
        early = cfg.build_cal1_settings(2, RANGES)     # built BEFORE the flag
        cfg.set_slurm_resource_env('1G', '2-00:00:00')
        for stage, build in builders().items():
            _, fixed, _ = build()
            check(f'{stage} followed the 2-day wall',
                  fixed['max_runtime'], dm.graceful_max_runtime('2-00:00:00'))
        _, fixed, _ = early()
        check('a builder made before the flag still follows it',
              fixed['max_runtime'], dm.graceful_max_runtime('2-00:00:00'))

        cfg.set_slurm_resource_env('1G', '00:30:00')
        _, fixed, _ = cfg.build_cal2_settings(2, RANGES)()
        check('and follows it down again',
              fixed['max_runtime'], dm.graceful_max_runtime('00:30:00'))
    finally:
        cfg.set_slurm_resource_env(os.environ['SIMPLICITY_SLURM_MEM'], before)
    check('the env is restored', os.environ['SIMPLICITY_SLURM_TIME'], before)


def test_it_survives_the_settings_round_trip():
    '''fixed_params is where a stated value wins over standard_values.json --
    the whole reason this is written rather than defaulted.'''
    print('\n=== and it is what a simulation actually reads')
    _, fixed, _ = cfg.build_cal2_settings(2, RANGES)()
    merged = {**sm.read_standard_parameters_values(), **fixed}
    check('the stated value wins over the standard values',
          merged['max_runtime'], fixed['max_runtime'])


def main():
    print(__doc__.strip().splitlines()[0])
    for section in (test_parse, test_margin, test_every_stage_states_it,
                    test_it_reaches_a_simulation, test_it_follows_the_flag,
                    test_it_survives_the_settings_round_trip):
        try:
            section()
        except Exception as exc:                                # noqa: BLE001
            print(f'  [FAIL] {section.__name__} raised '
                  f'{type(exc).__name__}: {exc}')
            failures.append(f'{section.__name__}: {type(exc).__name__}: {exc}')

    print('\n' + '=' * 70)
    if failures:
        print(f'{len(failures)} FAILED:')
        for line in failures:
            print(f'  {line}')
    else:
        print('the graceful stop tracks the walltime, and every stage states it')
    sys.exit(1 if failures else 0)


if __name__ == '__main__':
    main()
