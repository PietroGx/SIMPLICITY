#!/usr/bin/env python3
'''Do the paper figures read the output tree, and do they read ALL of it?

    python tests/test_figures_on_real_output.py --exp-num 2
    python tests/test_figures_on_real_output.py --exp-num 2 --exp-name impact_long_shedders_unbound
    python tests/test_figures_on_real_output.py --exp-num 2 --heavy --traceback

WHY

The v2.4.71 refactor moved every output directory a level deeper and renamed
all of them. nature_plots/ and long_paper_figures/ hold the largest block of
call sites it never touched, and they cannot be exercised against a synthetic
tree: they read PRODUCTION experiments by name, one per scenario.

WHY "IT RAN" IS NOT THE QUESTION

An earlier version of this file called each preprocessor and reported whether
it returned. Twelve did, and that told us almost nothing, because the thing
that would actually go wrong does not raise.

fig1.get_panel_a_data loops over scenarios with `except Exception: continue`.
A scenario whose read breaks is silently dropped and the function returns a
perfectly good DataFrame with one fewer curve in it. get_shared_valid_seeds
does the same per repeat. So a half-read experiment and a fully-read one look
identical from the outside -- which is this repo's standing failure mode, a
measurement quietly describing less than its name claims.

So every check here compares what came back against what the RECORD says
should be there -- settings.json, repeats.json, and the directories on disk --
and names the difference when they disagree.

HOW TO READ THE OUTPUT

  FAIL        the refactor's problem. A path helper, resolver or reader is
              wrong, or a preprocessor silently dropped something that exists.
  incomplete  it returned, but covered less than the record says exists. The
              missing items are named. This is the interesting column.
  no data     returned nothing, legitimately: an input that is absent, or a
              filter nothing passed. Counted apart.
  n/a         that experiment shape was not run at all.
'''
import argparse
import os
import sys
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'scripts', 'nature_plots'))
sys.path.insert(0, os.path.join(REPO, 'scripts', 'long_paper_figures'))

import matplotlib
matplotlib.use('Agg')

import simplicity.dir_manager as dm
import simplicity.jobs as jobs
import simplicity.output_manager as om
import simplicity.settings_manager as sm

failed, incomplete, nodata, notapplicable, fine = [], [], [], [], []


def describe(value):
    """What came back, by type -- not len() of whatever it is.

    The previous version printed `len(value) item(s)` for everything, so a
    returned PATH was reported as "139 item(s)": the length of the string.
    """
    import pandas as pd
    if value is None:
        return 'None'
    if isinstance(value, pd.DataFrame):
        return f'DataFrame {len(value)}x{len(value.columns)}'
    if isinstance(value, str):
        # only path-check things shaped like paths: a plain word is a status,
        # and reporting it as a MISSING path is noise
        if os.sep in value:
            return f'path {"exists" if os.path.exists(value) else "MISSING"}: ...{value[-52:]}'
        return value
    if isinstance(value, tuple):
        return f'tuple of {len(value)}: ' + ', '.join(describe(v) for v in value)
    if isinstance(value, dict):
        return f'dict, {len(value)} key(s)'
    if hasattr(value, '__len__'):
        return f'{type(value).__name__}, {len(value)} item(s)'
    return repr(value)[:60]


def probe(name, call, coverage=None):
    """Run one preprocessor, then ask whether it covered what exists.

    `coverage` returns (got, expected) as comparable sets, or None when there
    is nothing meaningful to compare against.
    """
    start = time.time()
    try:
        value = call()
    except Exception as exc:
        failed.append((name, f'{type(exc).__name__}: {exc}'))
        print(f'  [FAIL      ] {name}\n                {type(exc).__name__}: {exc}')
        if os.environ.get('FIGTEST_TRACEBACK'):
            traceback.print_exc()
        return None
    elapsed = time.time() - start
    summary = describe(value)

    empty = value is None or (hasattr(value, '__len__') and len(value) == 0)
    if empty:
        nodata.append((name, summary))
        print(f'  [no data   ] {name}: {summary}  ({elapsed:.1f}s)')
        return value

    if coverage is not None:
        try:
            got, expected = coverage(value)
        except Exception as exc:
            failed.append((name, f'coverage check raised {type(exc).__name__}: {exc}'))
            print(f'  [FAIL      ] {name}: coverage check raised {exc}')
            return value
        missing = set(expected) - set(got)
        extra = set(got) - set(expected)
        if missing or extra:
            detail = (f'covered {sorted(got)}; the record says {sorted(expected)}'
                      + (f'; MISSING {sorted(missing)}' if missing else '')
                      + (f'; UNEXPECTED {sorted(extra)}' if extra else ''))
            incomplete.append((name, detail))
            print(f'  [incomplete] {name}: {summary}\n                {detail}')
            return value
        print(f'  [ok        ] {name}: {summary}, covering all '
              f'{len(set(expected))}  ({elapsed:.1f}s)')
    else:
        print(f'  [ok        ] {name}: {summary}  ({elapsed:.1f}s)')
    fine.append(name)
    return value


# ----------------------------------------------------------- the record

def survey(exp_name, exp_num, scenarios):
    """What settings.json, repeats.json and the disk each say is there.

    Printed first because every check below is measured against it, and
    because a disagreement HERE is a refactor problem on its own -- the figure
    code has not even been reached yet.
    """
    print('\nthe record vs the disk')
    print(f'  {"scenario":<12}{"sims":>6}{"sods":>6}{"repeats":>9}'
          f'{"seed dirs":>11}  groups')
    present = {}
    for scenario in scenarios:
        experiment = f'{exp_name}_{scenario}_#{exp_num}'
        try:
            sims = sm.read_simulations(experiment)
            groups = [g['name'] for g in sm.read_groups(experiment)]
            sods = dm.get_simulation_output_dirs(experiment)
            ssods = [s for sod in sods
                     for s in dm.get_seeded_simulation_output_dirs(sod)]
            repeats = sum(len(jobs.read_repeats(experiment, g)) for g in groups)
        except Exception as exc:
            print(f'  {scenario:<12} -- {type(exc).__name__}: {exc}')
            continue
        print(f'  {scenario:<12}{len(sims):>6}{len(sods):>6}{repeats:>9}'
              f'{len(ssods):>11}  {groups}')
        if len(sims) != len(sods):
            failed.append((f'record vs disk [{scenario}]',
                           f'{len(sims)} simulation(s) in settings.json but '
                           f'{len(sods)} directories in 04_Output'))
        if repeats != len(ssods):
            failed.append((f'record vs disk [{scenario}]',
                           f'{repeats} repeat(s) in repeats.json but '
                           f'{len(ssods)} seed directories on disk'))
        present[scenario] = {'sods': sods, 'ssods': ssods, 'groups': groups,
                             'repeats': repeats}
    return present


def seed_dirs(present, scenario):
    return {int(dm.get_seed_from_SSOD(s))
            for s in present[scenario]['ssods']}


# ------------------------------------------------------------- figures

def figure_1(exp_num, exp_name, present):
    print('\nfigure 1')
    import fig1_preprocess_data as f1
    # each of these loops scenarios with `except Exception: continue`, so a
    # scenario that cannot be read vanishes instead of raising
    expected = {f1.get_clinical_label(s) for s in present}

    probe('fig1.get_panel_a_data',
          lambda: f1.get_panel_a_data(exp_num=exp_num, exp_name=exp_name),
          coverage=lambda df: (set(df['cohort']), expected))
    # this one DOES loop every scenario, and drops any that raises, so a
    # missing cohort is either a broken read or an empty input -- and those
    # look identical from the outside. Count the input per scenario so an
    # `incomplete` says which.
    panel_b = probe('fig1.get_panel_b_data',
                    lambda: f1.get_panel_b_data(exp_num=exp_num,
                                                exp_name=exp_name),
                    coverage=lambda df: (set(df['cohort']) if 'cohort' in df
                                         else set(), expected))
    if panel_b is not None and set(panel_b.get('cohort', [])) != expected:
        print('        why: usable rows per scenario (it needs individuals of '
              'the right type with a recorded end of infection)')
        for scenario in sorted(present):
            wanted = 'standard' if scenario == 'control' else 'long_shedder'
            usable = 0
            for ssod in present[scenario]['ssods']:
                try:
                    frame = om.read_individuals_data(ssod)
                    rows = frame[frame['type'] == wanted]
                    usable += int(rows['t_not_infected'].notna().sum())
                except Exception:
                    pass
            print(f'          {scenario:<10} {usable:>5} {wanted} row(s) with '
                  f't_not_infected')
    # control only, by contract: "standard individuals in control, divergence
    # from the outbreak root". Expecting every cohort here would be the test
    # being wrong, not the function.
    probe('fig1.get_model_global_clock',
          lambda: f1.get_model_global_clock(exp_num=exp_num, exp_name=exp_name),
          coverage=lambda df: (set(df['cohort']) if 'cohort' in df else set(),
                               {'Control'} if 'control' in present else set()))
    probe('fig1.get_model_intrahost_clock',
          lambda: f1.get_model_intrahost_clock(exp_num=exp_num,
                                               exp_name=exp_name))
    # these two read Data/RealWorldData -- tracked files, and the ones that
    # went missing in the quota cleanup
    probe('fig1.get_panel_c_data  (RealWorldData)', f1.get_panel_c_data)
    probe('fig1.get_panel_de_data (RealWorldData)', f1.get_panel_de_data)


def figure_2(exp_num, exp_name, present, heavy):
    print('\nfigure 2')
    import fig2_preprocess_data as f2
    scenarios = sorted(present)
    if 'control' not in present:
        print('  [n/a       ] fig2 selects its seeds from `control`, which has '
              'no output here')
        notapplicable.append(('fig2', 'no control experiment'))
        return
    control_seeds = seed_dirs(present, 'control')

    seeds = probe(
        'fig2.get_shared_valid_seeds',
        lambda: f2.get_shared_valid_seeds(exp_num, scenarios, exp_name=exp_name),
        # it samples from control's repeats that ran long enough, so it must be
        # a SUBSET of control's seeds -- a seed it returns that does not exist
        # on disk means it is reading the wrong tree
        coverage=lambda got: (set(got), set(got) & control_seeds))
    if not seeds:
        return
    print(f'        control has {len(control_seeds)} seed dir(s); '
          f'{len(seeds)} passed the run-length filter')

    seed, scenario = sorted(seeds)[0], scenarios[0]
    probe('fig2.get_target_seed_dir',
          lambda: f2.get_target_seed_dir(exp_num, scenario, seed,
                                         exp_name=exp_name),
          coverage=lambda p: ({os.path.basename(str(p))},
                              {f'seed_{seed:04d}'}))
    probe('fig2.get_fig2_freq_data',
          lambda: f2.get_fig2_freq_data(exp_num, scenario, seed,
                                        exp_name=exp_name))
    if heavy:
        probe('fig2.get_fig2_clustered_data',
              lambda: f2.get_fig2_clustered_data(exp_num, scenario, seed,
                                                 exp_name=exp_name))
        probe('fig2.get_fig2_divergence_data',
              lambda: f2.get_fig2_divergence_data(exp_num, scenario, [seed],
                                                  exp_name=exp_name))


def figure_3(exp_num, exp_name, present, heavy):
    print('\nfigure 3')
    import fig3_preprocess_data as f3
    scenario = sorted(present)[0]
    known = set(present[scenario]['sods'])
    probe('fig3._experiment_sod',
          lambda: f3._experiment_sod(exp_num, scenario, exp_name=exp_name),
          # it must hand back one of the simulation directories that exist,
          # not a path it built and never checked
          coverage=lambda p: ({str(p)}, {str(p)} & known))
    probe('fig3.get_panel_a_data', lambda: f3.get_panel_a_data(exp_num, scenario))
    if heavy:
        probe('fig3.get_panel_b_data',
              lambda: f3.get_panel_b_data(exp_num, sorted(present)))


def long_paper(exp_num):
    print('\nlong-paper figures  (a DIFFERENT experiment shape: the grid)')
    master = os.path.join(dm.get_data_dir(), f'master_grid_log_#{exp_num}.csv')
    if not os.path.isfile(master):
        print(f'  [n/a       ] no grid at {os.path.basename(master)} -- these '
              f'read the long-paper grid, not the impact pipeline')
        notapplicable.append(('long_paper', f'no {os.path.basename(master)}'))
        return
    import long_shedders_preprocess as lp
    probe('long_paper.read_master_log', lambda: lp.read_master_log(exp_num))
    probe('long_paper.get_baseline_sod', lambda: lp.get_baseline_sod(exp_num))


def library(exp_name, exp_num, present):
    print('\nplots_manager / output_manager, against one production experiment')
    import simplicity.plots_manager as pm
    scenario = sorted(present)[0]
    experiment = f'{exp_name}_{scenario}_#{exp_num}'
    sods = present[scenario]['sods']
    # no coverage check: this returns one row per intra-host LINEAGE, not per
    # simulation, so comparing its row indices against the simulation count
    # (which an earlier version did) compares two unrelated things and reports
    # a difference that means nothing.
    probe('om.get_IH_lineages_data_experiment',
          lambda: om.get_IH_lineages_data_experiment(experiment))
    # returns None; it is run for its side effect, so report that it did not
    # raise rather than letting describe() read the word "rendered" as a path
    probe('pm.plot_IH_lineage_distribution_grouped_by_simulation',
          lambda: (pm.plot_IH_lineage_distribution_grouped_by_simulation(
              experiment), 'no exception')[1])


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--exp-num', type=int, required=True)
    parser.add_argument('--exp-name', default='impact_long_shedders')
    parser.add_argument('--heavy', action='store_true',
                        help='also the SNP-matrix panels (minutes)')
    parser.add_argument('--traceback', action='store_true')
    args = parser.parse_args()
    if args.traceback:
        os.environ['FIGTEST_TRACEBACK'] = '1'

    import _scenarios
    scenarios = _scenarios.scenario_names(args.exp_name)
    print(f'scenarios in the config: {scenarios}')
    present = survey(args.exp_name, args.exp_num, scenarios)
    if not present:
        print('\nNo production output. Run the pipeline, or check --exp-num.')
        sys.exit(1)

    for stage in (lambda: figure_1(args.exp_num, args.exp_name, present),
                  lambda: figure_2(args.exp_num, args.exp_name, present, args.heavy),
                  lambda: figure_3(args.exp_num, args.exp_name, present, args.heavy),
                  lambda: long_paper(args.exp_num),
                  lambda: library(args.exp_name, args.exp_num, present)):
        try:
            stage()
        except Exception as exc:
            failed.append(('<stage setup>', f'{type(exc).__name__}: {exc}'))
            print(f'  [FAIL      ] the stage itself: {type(exc).__name__}: {exc}')
            if args.traceback:
                traceback.print_exc()

    print('\n' + '=' * 72)
    print(f'{len(fine)} ok, {len(incomplete)} incomplete, {len(nodata)} no data, '
          f'{len(notapplicable)} n/a, {len(failed)} FAILED')

    if nodata:
        print('\nno data (legitimate: an absent input, or a filter nothing passed)')
        for name, detail in nodata:
            print(f'  {name}: {detail}')
    if notapplicable:
        print('\nn/a (that experiment shape was not run)')
        for name, detail in notapplicable:
            print(f'  {name}: {detail}')
    if incomplete:
        print('\nINCOMPLETE -- returned, but covered less than exists:')
        for name, detail in incomplete:
            print(f'  {name}\n      {detail}')
    if failed:
        print('\nFAILED:')
        for name, detail in failed:
            print(f'  {name}\n      {detail}')
        print('\n--traceback for the full stack of each.')

    sys.exit(1 if (failed or incomplete) else 0)


if __name__ == '__main__':
    main()
