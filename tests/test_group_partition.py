#!/usr/bin/env python3
'''Stage 7c: a multi-group experiment is partitioned by the group each
simulation was ASSIGNED to, not by reading its parameters back and bucketing
on the float values.

    python tests/test_group_partition.py

WHY

Three readers listed a whole multi-group experiment and then rebuilt the
partition from parameter values:

  long_nsr_calibration_plot  bucketed on the (tau_3_long, R_long) floats, then
                             looked the group NAME up on those values rounded
                             to 3 decimals
  impact_long_shedders_cal_2 selected each scenario's rows by matching
                             (tau_3_long, long_shedders_ratio) against the
                             frozen table at 3 decimals
  nsr_sensitivity_analysis   rebuilt each (duration, k_v) cell from every
                             output directory, with groups that had no names
                             at all

All three produce a healthy fit either way. That is the point: a value-derived
partition is indistinguishable from a correct one until two groups happen to
collide, and then it silently fits their union and reports the result under one
of their names. CLAUDE.md's standing trap -- a measurement measuring something
other than what its name claims -- with the measurement being "this group".

WHAT IT ASSERTS

  1 equivalence   on the real cal_1 and cal_2 grids, the old value-derived
                  partition and the new group partition are the same partition.
                  A refactor that changes the numbers is a different change,
                  and this is what says it did not.
  2 divergence    on a grid where two groups share a cohort, the old way merges
                  them into one bucket and the new way refuses, by name.
  3 no leakage    each group's fit is fed only its own rows. The per-repeat
                  measurement is synthetic and group-dependent, so a fit that
                  saw a neighbour's rows lands on a different number -- which
                  real, noisy data could never show.
  4 single group  every production scenario builds exactly one group named
                  after itself. Sixteen call sites still list production
                  experiments unscoped; that is correct only while this holds.
  5 named cells   nsr_sensitivity's groups name themselves, one per cell, and
                  a group holding two cells raises.

Nothing is simulated: the settings, simulation-parameter and repeat files are
written by the real pipeline path, the output tree with them, and only the
per-repeat measurement is stubbed.
'''
import os
import shutil
import sys
import tempfile

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'scripts'))
sys.path.insert(0, os.path.join(REPO, 'scripts', 'experiments'))

import simplicity.dir_manager as dm
import simplicity.jobs as jobs
import simplicity.output_manager as om
import simplicity.runme as runme
import simplicity.settings_manager as sm
import simplicity.tuning.evolutionary_rate as er

failures = []


def check(label, got, want):
    ok = got == want
    print(f"  [{'ok ' if ok else 'FAIL'}] {label}: {got!r}")
    if not ok:
        failures.append(f'{label}\n        got  {got!r}\n        want {want!r}')
    return ok


def check_close(label, got, want, rtol=1e-3):
    ok = want != 0 and abs(got - want) / abs(want) < rtol
    print(f"  [{'ok ' if ok else 'FAIL'}] {label}: {got:.8g} (want {want:.8g})")
    if not ok:
        failures.append(f'{label}: got {got!r}, want {want!r} (rtol {rtol})')
    return ok


def check_raises(label, exc_type, call, *needles):
    try:
        call()
    except exc_type as exc:
        text = str(exc)
        missing = [n for n in needles if n not in text]
        if missing:
            print(f"  [FAIL] {label}: raised, but the message never names {missing}")
            print(f"         {text}")
            failures.append(f'{label}: message missing {missing}: {text}')
            return
        print(f"  [ok ] {label}: {exc_type.__name__} naming {list(needles)}")
        return
    except Exception as exc:                                    # noqa: BLE001
        print(f"  [FAIL] {label}: raised {type(exc).__name__}, not {exc_type.__name__}")
        print(f"         {exc}")
        failures.append(f'{label}: wrong exception {type(exc).__name__}: {exc}')
        return
    print(f"  [FAIL] {label}: did not raise")
    failures.append(f'{label}: did not raise {exc_type.__name__}')


# --------------------------------------------------------------------- setup
class NoSim:
    '''A runner that creates each repeat's output directory and simulates
    nothing. Everything upstream of it -- settings.json, 02_Simulations,
    03_Repeats, 04_Output/<group>/<stem>/parameters.json -- is written by the
    real path, so the readers under test walk the real tree.'''

    @staticmethod
    def run_seeded_simulations(experiment_name, _run_one, groups=None):
        for group, record in jobs.all_repeats(experiment_name):
            om.setup_output_directory(experiment_name, group, record)


class Scratch:
    '''An isolated Data/ for one section.'''

    def __init__(self):
        self.root = tempfile.mkdtemp(prefix='group_partition_')
        # both: sm._data_dir is a copy taken at import and does not follow
        # dm.set_data_dir
        dm.set_data_dir(self.root)
        sm._data_dir = self.root

    def build(self, name, settings_func):
        runme.run_experiment(name, settings_func, simplicity_runner=NoSim,
                             archive_experiment=False)
        return name

    def close(self):
        shutil.rmtree(self.root, ignore_errors=True)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def legacy_value_partition(experiment_name, parameters):
    '''The pre-7c bucketing, kept here as the reference the new code has to
    agree with: list the whole experiment, read these parameters back off each
    output directory, and group on the tuple of their values.'''
    buckets = {}
    for sod in dm.get_simulation_output_dirs(experiment_name):
        key = tuple(sm.get_parameter_value_from_simulation_output_dir(sod, p)
                    for p in parameters)
        buckets.setdefault(key, set()).add(os.path.basename(sod))
    return buckets


def group_partition(experiment_name):
    return {group: {os.path.basename(sod) for sod in
                    dm.get_simulation_output_dirs(experiment_name, group=group)}
            for group in dm.get_groups(experiment_name)}


def same_partition(left, right):
    '''Equal as partitions: same blocks, regardless of what they are keyed on.'''
    return (sorted(frozenset(block) for block in left.values()) ==
            sorted(frozenset(block) for block in right.values()))


# ------------------------------------------------------- synthetic measurement
#
# OSR = slope(group) * NSR, exactly, with a slope unique to each group. Fit a
# straight line to that and the inversion returns target/slope -- so the
# calibrated value a group reports IS the identity of the rows it was fed. One
# contaminating row from a neighbour moves it.
SLOPES = {}


def slope_for(group):
    return SLOPES.setdefault(group, 1.0 + 0.5 * len(SLOPES))


def install_stubs(nsr_parameter):
    '''Replace the per-repeat measurement with a deterministic function of the
    repeat's group. Returns a callable that puts the real ones back.'''
    real = {
        'read_final_time': om.read_final_time,
        'read_sequencing_data_regression': om.read_sequencing_data_regression,
        'extract_ih_regression_data': er.extract_ih_regression_data,
        'tempest_regression': er.tempest_regression,
    }

    def osr_of(ssod):
        sod = os.path.dirname(ssod)
        group = dm.get_group_from_SSOD(ssod)
        nsr = sm.get_parameter_value_from_simulation_output_dir(sod, nsr_parameter)
        return slope_for(group) * float(nsr)

    def frame(ssod, rows=40):
        return pd.DataFrame({'_osr': [osr_of(ssod)] * rows,
                             'individual_type': ['standard'] * rows})

    class Fit:
        def __init__(self, value):
            self.coef_ = [value]

    om.read_final_time = lambda ssod: 1000.0
    om.read_sequencing_data_regression = frame
    er.extract_ih_regression_data = frame
    er.tempest_regression = lambda df: Fit(float(df['_osr'].iloc[0]))

    def restore():
        om.read_final_time = real['read_final_time']
        om.read_sequencing_data_regression = real['read_sequencing_data_regression']
        er.extract_ih_regression_data = real['extract_ih_regression_data']
        er.tempest_regression = real['tempest_regression']

    return restore


def fit_result_files(experiment_name):
    directory = dm.get_experiment_fit_result_dir(experiment_name)
    return sorted(f for f in os.listdir(directory) if f.endswith('_fit_results.csv')
                  or '_fit_results_' in f)


# ===========================================================================
# 1 + 2 + 3  cal_1: one group per long-shedder cohort
# ===========================================================================
def section_cal_1():
    print('\n=== cal_1: partition, divergence, and what each fit was fed')
    import impact_long_shedders_unbound_config as cfg
    import long_nsr_calibration_plot as plot

    with Scratch() as scratch:
        name = scratch.build(
            'cal1_probe',
            cfg.build_cal1_settings(2, {'min': 1e-4, 'max': 1e-2, 'steps': 3}))

        groups = plot._cal1_groups(name)
        check('cal_1 groups, in declared order', [g['name'] for g in groups],
              ['SOT', 'HIV_low+HIV_high', 'edge_case'])

        print('\n  the old partition and the new one are the same partition')
        legacy = legacy_value_partition(name, ('tau_3_long', 'R_long'))
        current = group_partition(name)
        check('blocks: value-derived vs assigned', len(legacy), len(current))
        check('and they hold the same simulations',
              same_partition(legacy, current), True)

        print('\n  each cohort is one (tau_3_long, R_long)')
        check('cohorts', sorted((g['tau_3_long'], g['R_long']) for g in groups),
              sorted(cfg.unique_long_tau_r_long_pairs(
                  sm.read_standard_parameters_values())))

        print('\n  each fit saw only its own group')
        restore = install_stubs('nucleotide_substitution_rate_long')
        try:
            results = plot.plot_and_fit_long_nsr_calibration(
                name, target_osr_long=0.002, model_type='lin',
                min_seq=10, min_len=10)
        finally:
            restore()

        check('one calibrated value per group', len(results), len(groups))
        for group in groups:
            key = (group['tau_3_long'], group['R_long'])
            check_close(f"  {group['name']} inverts to target/slope",
                        results[key], 0.002 / slope_for(group['name']))

        print('\n  and the saved artefacts carry the recorded names')
        saved = fit_result_files(name)
        check('one fit-results file per group', len(saved), len(groups))
        for group in groups:
            check(f"  {group['name']} has its own file",
                  any(f.endswith(f"_{group['name']}.csv") for f in saved), True)

        frozen = pd.read_csv(plot._calibrated_long_nsr_path(name))
        check('the frozen table names its groups',
              sorted(frozen['group'].tolist()),
              sorted(g['name'] for g in groups))

    print('\n  two groups on one cohort: merged before, refused now')
    with Scratch() as scratch:
        # Same cohort, different k_v -- the exact shape nsr_sensitivity sweeps.
        # Bucketing on (tau_3_long, R_long) cannot see the difference.
        shared = {'tau_3_long': 100.0, 'R_long': 1.1,
                  'nucleotide_substitution_rate_long': [1e-4, 1e-3]}
        name = scratch.build('cal1_collision', lambda: (
            {'_scenario_groups': [
                dict(shared, name='cohort_kv_low', IH_virus_emergence_rate=0.01),
                dict(shared, name='cohort_kv_high', IH_virus_emergence_rate=0.1)]},
            {'population_size': 20, 'final_time': 2}, 1))

        legacy = legacy_value_partition(name, ('tau_3_long', 'R_long'))
        current = group_partition(name)
        check('the old way sees one block where there are two',
              (len(legacy), len(current)), (1, 2))
        check('so the two partitions disagree',
              same_partition(legacy, current), False)
        check_raises('and _cal1_groups refuses instead of picking one',
                     RuntimeError, lambda: plot._cal1_groups(name),
                     'cohort_kv_low', 'cohort_kv_high')


# ===========================================================================
# 1 + 3  cal_2: one group per scenario
# ===========================================================================
def section_cal_2():
    print('\n=== cal_2: one group per scenario, matched by name')
    import impact_long_shedders_config as cfg
    import impact_long_shedders_cal_2 as cal_2

    sp = sm.read_standard_parameters_values()
    long_nsr_by_group = {pair: 5e-4
                         for pair in cfg.unique_long_tau_r_long_pairs(sp)}
    ranges = {s['name']: {'min': 1e-4, 'max': 1e-2, 'steps': 3}
              for s in cfg.SCENARIOS}

    with Scratch() as scratch:
        groups, frozen_rows = cfg.build_cal2_scenario_groups(
            ranges, long_nsr_by_group)
        name = scratch.build('cal2_probe',
                             cfg.build_cal2_settings(groups, 2))

        expected = [row['scenario_name'] for row in frozen_rows]
        check('groups are the scenarios', dm.get_groups(name), expected)

        print('\n  the old partition and the new one are the same partition')
        legacy = legacy_value_partition(
            name, ('tau_3_long', 'long_shedders_ratio'))
        current = group_partition(name)
        check('blocks: value-derived vs assigned', len(legacy), len(current))
        check('and they hold the same simulations',
              same_partition(legacy, current), True)

        print('\n  each scenario fit saw only its own group')
        SLOPES.clear()
        restore = install_stubs('nucleotide_substitution_rate')
        try:
            calibrated = cal_2.compute_standard_nsr_per_scenario(
                name, frozen_rows, target_osr_std=0.002, model_type='lin',
                min_seq=10, min_len=10)
        finally:
            restore()

        check('one calibrated value per scenario',
              sorted(calibrated), sorted(expected))
        for scenario in expected:
            check_close(f'  {scenario} inverts to target/slope',
                        calibrated[scenario], 0.002 / slope_for(scenario))

        saved = fit_result_files(name)
        check('one fit-results file per scenario', len(saved), len(expected))

        print('\n  a frozen table the run does not match is named, not guessed')
        stranger = list(frozen_rows) + [dict(frozen_rows[0],
                                             scenario_name='not_run')]
        check_raises('missing group reported by name', RuntimeError,
                     lambda: cal_2.compute_standard_nsr_per_scenario(
                         name, stranger, 0.002, 'lin', 10, 10),
                     'not_run')


# ===========================================================================
# 4  production is one group, which is why 16 sites may list it unscoped
# ===========================================================================
def section_production_is_single_group():
    print('\n=== production: one experiment, one group, named after the scenario')
    import impact_long_shedders_unbound_config as unbound
    import impact_long_shedders_config as bound

    row = {
        'scenario_name': 'HIV_high', 'tau_3_long': 100.0,
        'long_shedders_ratio': 0.12, 'susceptibility_long': 1.0, 'R': 1.05,
        'IH_virus_emergence_rate': 0.1, 'R_long': 1.1,
        'nucleotide_substitution_rate_long': 5e-4,
        'nucleotide_substitution_rate': 1.4e-4,
    }
    control = dict(row, scenario_name='control', long_shedders_ratio=0.0,
                   nucleotide_substitution_rate_long=np.nan)

    for label, build in (('unbound', unbound.build_exp_scenario_settings),
                         ('bound', bound.build_exp_scenario_settings)):
        for case in (row, control):
            varying, _fixed, _seeds = build(pd.Series(case), 5)()
            declared = varying['_scenario_groups']
            check(f"{label}/{case['scenario_name']}: exactly one group",
                  len(declared), 1)
            check(f"{label}/{case['scenario_name']}: named after the scenario",
                  declared[0]['name'], case['scenario_name'])

    with Scratch() as scratch:
        name = scratch.build(
            'impact_long_shedders_unbound_HIV_high_probe',
            unbound.build_exp_scenario_settings(pd.Series(row), 2))
        check('and on disk it is one group', dm.get_groups(name), ['HIV_high'])
        check('so an unscoped listing is that group',
              dm.get_simulation_output_dirs(name),
              dm.get_simulation_output_dirs(name, group='HIV_high'))


# ===========================================================================
# 5  nsr_sensitivity: cells name themselves
# ===========================================================================
def section_nsr_sensitivity_cells():
    print('\n=== nsr_sensitivity: one named group per (duration, k_v) cell')
    import nsr_sensitivity_analysis as nsr

    names = [nsr.cell_group_name(d, k)
             for d in nsr.DURATIONS for k in nsr.K_V_VALUES]
    check('one name per cell, all distinct', len(set(names)), len(names))
    check('and none of them could become a path',
          [n for n in names if '/' in n or os.sep in n], [])

    with Scratch() as scratch:
        sp = sm.read_standard_parameters_values()
        durations, k_vs = nsr.QUICK_DURATIONS, nsr.QUICK_K_V
        name = scratch.build('nsr_cells_probe', nsr.build_settings(
            durations, k_vs, [1e-4, 1e-3], 1, 40, 4, sp))

        cells = nsr.cell_groups(name, sp)
        check('a group per cell', len(cells), len(durations) * len(k_vs))
        check('named as the grid names them', sorted(cells),
              sorted(nsr.cell_group_name(d, k) for d in durations for k in k_vs))
        check('each cell reports its own duration and k_v',
              sorted((c['inf_duration'], c['k_v']) for c in cells.values()),
              sorted((d, k) for d in durations for k in k_vs))

        print('\n  the old partition and the new one are the same partition')
        legacy = legacy_value_partition(
            name, ('tau_3_long', 'IH_virus_emergence_rate'))
        check('blocks: value-derived vs assigned',
              same_partition(legacy, group_partition(name)), True)

    print('\n  a group holding two cells is refused')
    with Scratch() as scratch:
        name = scratch.build('nsr_cells_mixed', lambda: (
            {'_scenario_groups': [{'name': 'two_cells',
                                   'tau_3_long': [50.0, 90.0],
                                   'IH_virus_emergence_rate': 0.1,
                                   'final_time': 365}]},
            {'population_size': 20}, 1))
        check_raises('cell_groups names the parameter that varies',
                     RuntimeError,
                     lambda: nsr.cell_groups(name,
                                             sm.read_standard_parameters_values()),
                     'two_cells', 'tau_3_long')


# ===========================================================================
# 6  the cluster checker itself: can it tell "unchanged" from "changed"?
# ===========================================================================
def section_calibration_checker():
    """check_calibration_unchanged.py is the one thing here that runs against a
    REAL finished run, on the cluster, where this suite cannot follow it. So
    exercise its comparison on a synthetic run instead: a checker that cannot
    detect a change is worse than no checker, and that is exactly what a clean
    'reproduces exactly' from an untested script would look like."""
    print('\n=== the cluster checker detects a change, and only a change')
    import check_calibration_unchanged as checker
    import impact_long_shedders_unbound_config as cfg
    import long_nsr_calibration_plot as plot

    with Scratch() as scratch:
        name = scratch.build(
            'impact_long_shedders_unbound_cal_long_#77',
            cfg.build_cal1_settings(2, {'min': 1e-4, 'max': 1e-2, 'steps': 3}))

        SLOPES.clear()
        restore = install_stubs('nucleotide_substitution_rate_long')
        try:
            # a finished run: the fit, the frozen CSV, the fit-results files
            plot.plot_and_fit_long_nsr_calibration(
                name, 0.002, model_type='lin', min_seq=10, min_len=10)
            saved = sorted(os.listdir(dm.get_experiment_fit_result_dir(name)))

            print('\n  re-fitting the same output reproduces it')
            checker.failures.clear()
            checker.skipped.clear()
            checker.stage_1(
                {'cal_1': 'impact_long_shedders_unbound_cal_long'}, 77,
                rtol=1e-9, keep=False)
            check('no differences reported', checker.failures, [])
            check('and nothing was skipped', checker.skipped, [])
            check('the run\'s own Fit_results are back',
                  sorted(os.listdir(dm.get_experiment_fit_result_dir(name))),
                  saved)

            print('\n  and a changed measurement is caught')
            for group in list(SLOPES):
                SLOPES[group] *= 1.5          # every cohort now calibrates 1/1.5x
            checker.failures.clear()
            checker.skipped.clear()
            checker.stage_1(
                {'cal_1': 'impact_long_shedders_unbound_cal_long'}, 77,
                rtol=1e-9, keep=False)
            check('one difference reported per cohort',
                  len(checker.failures), len(SLOPES))
            checker.failures.clear()
            checker.skipped.clear()
        finally:
            restore()


def main():
    section_cal_1()
    section_cal_2()
    section_production_is_single_group()
    section_nsr_sensitivity_cells()
    section_calibration_checker()

    print('\n' + '=' * 70)
    if failures:
        print(f'{len(failures)} FAILURE(S):')
        for line in failures:
            print(f'  {line}')
    else:
        print('every multi-group reader partitions by the assigned group, and '
              'agrees with\nthe old value-derived partition wherever that one '
              'was still right')
    sys.exit(1 if failures else 0)


if __name__ == '__main__':
    main()
