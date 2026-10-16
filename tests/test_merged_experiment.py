#!/usr/bin/env python3
'''Stage 8e-8g: a run is ONE experiment whose phases are named groups.

    python tests/test_merged_experiment.py

WHY

Before 8e a run wrote nine experiment directories -- one per calibration
stage, one per scenario, plus setup data and sanity. Each was homogeneous, so
listing one without naming a group could not reach anything else, and every
reader did exactly that.

Now the long calibration's cohorts, the standard sweep and all five scenarios
live in impact_long_shedders_#N together. An unscoped listing returns the lot.
Every one of those simulations has a tau_3_long, an R_long and a
nucleotide_substitution_rate, so a reader that forgets to scope does not
crash -- it fits a calibration sweep and reports the number under a scenario's
name. That is CLAUDE.md's standing trap, and the figures have been bitten by
its shape twice already.

So this does not check that the phases "work". It checks the two things that
make the merge safe:

  1 STRUCTURE   three phases build one experiment of nine groups; each phase
                runs only its own; ids continue across phases so no two
                simulations share an output directory; each group keeps its
                own n_seeds.
  2 SCOPING     every reader that takes an experiment name returns rows from
                the group it claims and from nothing else -- checked by
                giving the calibration groups a value no scenario has and
                requiring it never to appear in a scenario's rows.

Nothing is simulated: the runner is stubbed and only the output tree is made.
'''
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
sys.path.insert(0, os.path.join(REPO, 'scripts', 'nature_plots'))

import simplicity.dir_manager as dm
import simplicity.jobs as jobs
import simplicity.output_manager as om
import simplicity.runme as runme
import simplicity.settings_manager as sm

import impact_long_shedders_config as cfg

EXP_NUM = 1
NAME = f'{cfg.EXP_NAME}_#{EXP_NUM}'
CAL_SEEDS, EXP_SEEDS = 2, 3
RANGES = {'min': 1e-4, 'max': 1e-2, 'steps': 3}

failures = []


def check(label, got, want):
    ok = got == want
    print(f"  [{'ok ' if ok else 'FAIL'}] {label}: {got!r}")
    if not ok:
        failures.append(f'{label}: got {got!r}, want {want!r}')


class Runner:
    '''Builds the output tree a real run would leave, and records the scope
    it was given.'''
    calls = []

    @staticmethod
    def run_seeded_simulations(experiment_name, run_one, groups=None):
        Runner.calls.append(list(groups) if groups else None)
        for group, record in jobs.all_repeats(experiment_name, groups):
            om.setup_output_directory(experiment_name, group, record)


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


def build():
    '''The three phases, in the order and with the append flags the real
    stage scripts use.'''
    runme.run_experiment_phase(NAME, cfg.build_cal1_settings(CAL_SEEDS, RANGES),
                               Runner, append=False)
    runme.run_experiment_phase(NAME, cfg.build_cal2_settings(CAL_SEEDS, RANGES),
                               Runner, append=True)
    runme.run_experiment_phase(NAME, cfg.build_exp_phase_settings(frozen_rows(),
                                                                  EXP_SEEDS),
                               Runner, append=True)


# ===========================================================================
def test_structure():
    print('\n=== one experiment, nine groups')
    groups = dm.get_groups(NAME)
    expected = (
        [cfg.calib_long_group(c)
         for c in ('SOT', 'HIV_low+HIV_high', 'edge_case')]
        + [cfg.calib_std_group()]
        + [cfg.exp_group(s['name']) for s in cfg.SCENARIOS])
    check('every phase is a group of one experiment, in order',
          groups, expected)

    check('each phase ran only its own groups', Runner.calls,
          [expected[:3], expected[3:4], expected[4:]])

    seeds = {g: sm.read_group_n_seeds(NAME, g) for g in groups}
    check('calibration groups keep --cal-seeds',
          sorted({seeds[g] for g in expected[:4]}), [CAL_SEEDS])
    check('experiment groups keep --exp-seeds',
          sorted({seeds[g] for g in expected[4:]}), [EXP_SEEDS])

    sims = sm.read_simulations(NAME)
    ids = [record['id'] for record in sims]
    check('simulation ids continue across phases', ids, list(range(len(ids))))
    standard = sm.read_standard_parameters_values()
    stems = {sm.simulation_stem(record, standard) for record in sims}
    check('every output directory is distinct -- nothing is overwritten',
          len(stems), len(sims))


# ===========================================================================
def test_an_unscoped_listing_really_does_span_phases():
    '''The premise. If this ever stops holding, the scoping below is moot and
    the checks would pass for the wrong reason.'''
    print('\n=== an unscoped listing spans every phase')
    everything = dm.get_simulation_output_dirs(NAME)
    one_group = dm.get_simulation_output_dirs(NAME, group=cfg.exp_group('SOT'))
    check('unscoped returns more than one group holds',
          len(everything) > len(one_group), True)
    check('and it reaches the calibration sweep',
          any(cfg.CALIB_LONG_PREFIX in path.split(os.sep)[-3:][0]
              or cfg.CALIB_LONG_PREFIX in path for path in everything), True)


# ===========================================================================
def test_figures_read_only_their_scenario():
    '''The 8f guard, on the resolver every panel goes through.'''
    print('\n=== a scenario resolves to its own group and nothing else')
    for scenario in (s['name'] for s in cfg.SCENARIOS):
        group = cfg.exp_group(scenario)
        sods = cfg.scenario_output_dirs(cfg.EXP_NAME, scenario, EXP_NUM)
        check(f'{scenario}: every directory is in its group',
              all(group in sod.split(os.sep) for sod in sods), True)
        check(f'{scenario}: one grid point, as the frozen table has one row',
              len(sods), 1)

    # the calibration cohorts share SOT's and edge_case's names after the
    # prefix; without it they would be the same group
    sot = cfg.scenario_output_dirs(cfg.EXP_NAME, 'SOT', EXP_NUM)
    calib = dm.get_simulation_output_dirs(NAME,
                                          group=cfg.calib_long_group('SOT'))
    check('the SOT scenario and the SOT cohort are different directories',
          set(sot) & set(calib), set())
    check('and the cohort is the one with the swept grid', len(calib), 3)

    check('an unknown scenario resolves to nothing, rather than everything',
          cfg.scenario_output_dirs(cfg.EXP_NAME, 'not_a_scenario', EXP_NUM), [])
    check('has_output agrees',
          [cfg.scenario_output_dirs(cfg.EXP_NAME, 'SOT', EXP_NUM) != [],
           cfg.scenario_output_dirs(cfg.EXP_NAME, 'SOT', 999) != []],
          [True, False])


# ===========================================================================
def test_readers_that_take_a_whole_experiment():
    '''8g: the two deliberate all-group listings, and the n_seeds read that
    stage 8a made raise when groups disagree.'''
    print('\n=== whole-experiment readers')

    import check_completed_simulations as checker

    # read_n_seeds_file is what this used to call; the merged experiment has
    # two different seed counts, so it must refuse rather than pick one
    try:
        sm.read_n_seeds_file(NAME)
        check('read_n_seeds_file refuses a mixed experiment', False, True)
    except Exception as exc:                                    # noqa: BLE001
        named = (cfg.calib_std_group() in str(exc)
                 and cfg.exp_group('SOT') in str(exc))
        check('read_n_seeds_file refuses a mixed experiment, naming groups',
              named, True)

    # ... so the checker has to go group by group, against each own n_seeds
    import contextlib
    import io as _io
    buf = _io.StringIO()
    with contextlib.redirect_stdout(buf):
        checker.count_and_analyze_simulations(NAME)
    out = buf.getvalue()
    check('the completeness check reports every group',
          all(f'group {g}' in out for g in dm.get_groups(NAME)), True)
    check('and against the right n_seeds per phase',
          (f'n_seeds={CAL_SEEDS}' in out and f'n_seeds={EXP_SEEDS}' in out),
          True)

    buf = _io.StringIO()
    with contextlib.redirect_stdout(buf):
        checker.count_and_analyze_simulations(NAME, [cfg.exp_group('SOT')])
    only = buf.getvalue()
    check('--group narrows it to one phase',
          (f'group {cfg.exp_group("SOT")}' in only
           and cfg.calib_std_group() not in only), True)


# ===========================================================================
def test_cal2_fit_names_the_group_it_was_scoped_to():
    """The stage that actually broke on the HPC.

    compute_standard_nsr reads its rows scoped, then has to NAME the fit for
    the results filename. It took that name from dm.get_groups(numbered) and
    required exactly one -- true while cal_2 owned an experiment, false since
    8e, where cal_1's cohorts are in the same experiment. So the stage died
    AFTER its 300 simulations had run, on an experiment that was fine.

    Only the per-seed readers and the regressor are stubbed; the scoping, the
    group listing and the naming are the real ones.
    """
    print('\n=== cal_2 fits under its own group, not whatever is on disk')
    import impact_long_shedders_cal_2 as c2

    class Fit:
        def eval(self, x=None):
            return [1.0 for _ in x]

    seen = {}
    saved = (c2.om.read_final_time, c2.om.read_sequencing_data_regression,
             c2.er.tempest_regression,
             c2.er.fit_observed_substitution_rate_regressor,
             c2.er.compute_calibrated_parameter)

    def record_fit(numbered, df, model, parameter_name=None,
                   experiment_group=None):
        seen['group'] = experiment_group
        seen['rows'] = len(df)
        return Fit()

    c2.om.read_final_time = lambda ssod: 365.0
    c2.om.read_sequencing_data_regression = lambda ssod: pd.DataFrame(
        {'individual_type': ['standard'] * 40})
    c2.er.tempest_regression = lambda df: type('R', (), {'coef_': [0.0013]})()
    c2.er.fit_observed_substitution_rate_regressor = record_fit
    c2.er.compute_calibrated_parameter = lambda *a, **k: 7e-5
    try:
        nsr = c2.compute_standard_nsr(NAME, 0.0013)
        check('it runs at all on a merged experiment', isinstance(nsr, float),
              True)
        check('and names the calib_std group, not one of the four on disk',
              seen.get('group'), cfg.calib_std_group())
        # one row per SEED: calib_std's 3 grid points x CAL_SEEDS. Unscoped
        # it would pull every seed of all 17 simulations in the experiment.
        check('the rows it fitted came only from that group',
              seen.get('rows'), 3 * CAL_SEEDS)

        # and an explicit scope still wins
        seen.clear()
        c2.compute_standard_nsr(NAME, 0.0013, groups=[cfg.calib_std_group()])
        check('an explicit scope names the same group',
              seen.get('group'), cfg.calib_std_group())

        # the guard that remains: a fit is one group, and it says so
        try:
            c2.compute_standard_nsr(NAME, 0.0013, groups=[
                cfg.calib_std_group(), cfg.calib_long_group('SOT')])
            check('two groups for one fit is refused', False, True)
        except RuntimeError as exc:
            check('two groups for one fit is refused, naming them',
                  'one fit needs one group' in str(exc), True)
    finally:
        (c2.om.read_final_time, c2.om.read_sequencing_data_regression,
         c2.er.tempest_regression,
         c2.er.fit_observed_substitution_rate_regressor,
         c2.er.compute_calibrated_parameter) = saved


def main():
    print(__doc__.strip().splitlines()[0])
    root = tempfile.mkdtemp(prefix='merged_exp_')
    data = os.path.join(root, 'Data')
    os.makedirs(data)
    saved = os.getcwd()
    os.chdir(root)
    dm.set_data_dir(data)
    sm._data_dir = data
    try:
        build()
        for section in (test_structure,
                        test_an_unscoped_listing_really_does_span_phases,
                        test_figures_read_only_their_scenario,
                        test_readers_that_take_a_whole_experiment,
                        test_cal2_fit_names_the_group_it_was_scoped_to):
            # a section that RAISES is a failure, not a crash: the scoping
            # guards raise by design, and a traceback here would end the run
            # with the later sections unreported
            try:
                section()
            except Exception as exc:                            # noqa: BLE001
                print(f'  [FAIL] {section.__name__} raised '
                      f'{type(exc).__name__}: {exc}')
                failures.append(f'{section.__name__}: '
                                f'{type(exc).__name__}: {exc}')
    finally:
        os.chdir(saved)
        shutil.rmtree(root, ignore_errors=True)

    print('\n' + '=' * 70)
    if failures:
        print(f'{len(failures)} FAILED:')
        for line in failures:
            print(f'  {line}')
    else:
        print('one experiment, phases as groups, and every reader stays '
              'inside the group it names')
    sys.exit(1 if failures else 0)


if __name__ == '__main__':
    main()
