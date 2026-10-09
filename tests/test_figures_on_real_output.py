#!/usr/bin/env python3
'''Do the paper figures still read the output tree? Run this ON THE CLUSTER.

    python tests/test_figures_on_real_output.py --exp-num 2
    python tests/test_figures_on_real_output.py --exp-num 2 \\
        --exp-name impact_long_shedders_unbound

WHY

The v2.4.71 data flow refactor moved every output directory a level deeper and
renamed all of them. The calibration and sanity figures are verified -- the
cluster pipeline produced them. The PAPER figures are not: nature_plots/ and
long_paper_figures/ have never run against the new layout, and between them
they hold the largest block of untouched get_simulation_output_dirs and
get_parameter_value_from_simulation_output_dir call sites in the repo.

They cannot be exercised against a synthetic tree: they read PRODUCTION
experiments by name, one per scenario (impact_long_shedders_<scenario>_#<n>),
and expect real sequence and lineage data. So this runs against whatever is
actually on disk.

WHAT IT DOES

Calls the real preprocessors -- the functions the figure scripts call -- and
reports, per entry point, whether it returned, raised, or had nothing to read.
It renders nothing and writes nothing: a preprocessor that returns its data is
the thing in question, and plotting adds minutes without adding an answer.

IT IS A REPORT, NOT A PASS/FAIL SUITE. "no data" is a legitimate outcome for a
small or partial run and is counted separately from a failure. What matters is
the FAILED column: an exception from a path helper, a resolver or a reader is
the refactor's problem. An exception about too few points, an empty frame or a
missing optional input is the data's.

Exit code is non-zero only if something actually raised.
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
sys.path.insert(0, os.path.join(REPO, 'scripts', 'experiments'))

import matplotlib
matplotlib.use('Agg')

import simplicity.dir_manager as dm

results = []          # (name, outcome, detail, seconds)


def run(name, call):
    """Call one preprocessor and classify what came back."""
    start = time.time()
    try:
        value = call()
    except Exception as exc:
        elapsed = time.time() - start
        results.append((name, 'FAILED', f'{type(exc).__name__}: {exc}', elapsed))
        print(f'  [FAILED ] {name}\n             {type(exc).__name__}: {exc}')
        if os.environ.get('FIGTEST_TRACEBACK'):
            traceback.print_exc()
        return None
    elapsed = time.time() - start
    size = len(value) if hasattr(value, '__len__') else None
    if value is None or size == 0:
        results.append((name, 'no data', f'returned {value!r:.40}', elapsed))
        print(f'  [no data] {name}  ({elapsed:.1f}s)')
    else:
        detail = f'{size} item(s)' if size is not None else repr(value)[:60]
        results.append((name, 'ok', detail, elapsed))
        print(f'  [ok     ] {name}: {detail}  ({elapsed:.1f}s)')
    return value


def which_experiments(exp_name, exp_num):
    """What production output actually exists, via the scripts' own helper."""
    import _scenarios
    print(f'\nscenario discovery  ({exp_name}_#{exp_num})')
    try:
        names = _scenarios.scenario_names(exp_name)
        print(f'  scenarios in the config : {names}')
    except Exception as exc:
        print(f'  [FAILED ] scenario_names: {type(exc).__name__}: {exc}')
        return []
    present = []
    for scenario in names:
        experiment = f'{exp_name}_{scenario}_#{exp_num}'
        try:
            sods = dm.get_simulation_output_dirs(experiment)
            groups = dm.get_groups(experiment)
            print(f'  {scenario:<10} {len(sods)} simulation(s), group(s) {groups}')
            present.append(scenario)
        except Exception as exc:
            print(f'  {scenario:<10} -- {type(exc).__name__}: {exc}')
    return present


def figure_1(exp_num, exp_name, scenarios):
    print('\nfigure 1 preprocessors')
    import fig1_preprocess_data as f1
    run('fig1.get_panel_a_data',
        lambda: f1.get_panel_a_data(exp_num=exp_num, exp_name=exp_name))
    run('fig1.get_panel_b_data',
        lambda: f1.get_panel_b_data(exp_num=exp_num, exp_name=exp_name))
    # panels c and d/e read Data/RealWorldData -- tracked files, and the ones
    # that went missing in the quota cleanup
    run('fig1.get_panel_c_data  (RealWorldData)', f1.get_panel_c_data)
    run('fig1.get_panel_de_data (RealWorldData)', f1.get_panel_de_data)
    run('fig1.get_model_global_clock',
        lambda: f1.get_model_global_clock(exp_num=exp_num, exp_name=exp_name))
    run('fig1.get_model_intrahost_clock',
        lambda: f1.get_model_intrahost_clock(exp_num=exp_num, exp_name=exp_name))


def figure_2(exp_num, exp_name, scenarios, heavy):
    print('\nfigure 2 preprocessors')
    import fig2_preprocess_data as f2
    seeds = run('fig2.get_shared_valid_seeds',
                lambda: f2.get_shared_valid_seeds(exp_num, scenarios,
                                                  exp_name=exp_name))
    if not seeds:
        print('        (no shared valid seed -- the per-seed panels need one)')
        return
    seed = sorted(seeds)[0]
    scenario = scenarios[0]
    print(f'        using scenario={scenario} seed={seed}')
    run('fig2.get_target_seed_dir',
        lambda: f2.get_target_seed_dir(exp_num, scenario, seed,
                                       exp_name=exp_name))
    run('fig2.get_fig2_freq_data',
        lambda: f2.get_fig2_freq_data(exp_num, scenario, seed,
                                      exp_name=exp_name))
    if heavy:
        run('fig2.get_fig2_clustered_data',
            lambda: f2.get_fig2_clustered_data(exp_num, scenario, seed,
                                               exp_name=exp_name))
        run('fig2.get_fig2_divergence_data',
            lambda: f2.get_fig2_divergence_data(exp_num, scenario, [seed],
                                                exp_name=exp_name))


def figure_3(exp_num, exp_name, scenarios, heavy):
    print('\nfigure 3 preprocessors')
    import fig3_preprocess_data as f3
    group = scenarios[0]
    run('fig3._experiment_sod',
        lambda: f3._experiment_sod(exp_num, group, exp_name=exp_name))
    run('fig3.get_panel_a_data',
        lambda: f3.get_panel_a_data(exp_num, group))
    if heavy:
        # panel b and c build SNP matrices; slow, and the pair that compared
        # 1,494 long-shedder lineages against 6 sequenced standards
        run('fig3.get_panel_b_data',
            lambda: f3.get_panel_b_data(exp_num, scenarios))


def long_paper(exp_num):
    print('\nlong-paper preprocessors  (a different experiment shape)')
    import long_shedders_preprocess as lp
    run('long_paper.read_master_log', lambda: lp.read_master_log(exp_num))
    run('long_paper.get_baseline_sod', lambda: lp.get_baseline_sod(exp_num))


def library(exp_name, exp_num, scenarios):
    print('\nplots_manager entry points, against one production experiment')
    import simplicity.plots_manager as pm
    import simplicity.output_manager as om
    experiment = f'{exp_name}_{scenarios[0]}_#{exp_num}'
    run('om.get_IH_lineages_data_experiment',
        lambda: om.get_IH_lineages_data_experiment(experiment))
    run('pm.plot_IH_lineage_distribution_grouped_by_simulation',
        lambda: pm.plot_IH_lineage_distribution_grouped_by_simulation(experiment)
        or 'rendered')


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--exp-num', type=int, required=True)
    parser.add_argument('--exp-name', default='impact_long_shedders',
                        help='impact_long_shedders or impact_long_shedders_unbound')
    parser.add_argument('--heavy', action='store_true',
                        help='also run the SNP-matrix panels (minutes)')
    parser.add_argument('--traceback', action='store_true')
    args = parser.parse_args()
    if args.traceback:
        os.environ['FIGTEST_TRACEBACK'] = '1'

    scenarios = which_experiments(args.exp_name, args.exp_num)
    if not scenarios:
        print('\nNo production output found. Nothing to test against -- run the '
              'pipeline first, or check --exp-num / --exp-name.')
        sys.exit(1)

    for stage in (lambda: figure_1(args.exp_num, args.exp_name, scenarios),
                  lambda: figure_2(args.exp_num, args.exp_name, scenarios, args.heavy),
                  lambda: figure_3(args.exp_num, args.exp_name, scenarios, args.heavy),
                  lambda: long_paper(args.exp_num),
                  lambda: library(args.exp_name, args.exp_num, scenarios)):
        try:
            stage()
        except Exception as exc:
            print(f'  [FAILED ] the stage itself: {type(exc).__name__}: {exc}')
            results.append(('<stage setup>', 'FAILED',
                            f'{type(exc).__name__}: {exc}', 0.0))

    failed = [r for r in results if r[1] == 'FAILED']
    nodata = [r for r in results if r[1] == 'no data']
    ok = [r for r in results if r[1] == 'ok']

    print('\n' + '=' * 72)
    print(f'{len(ok)} ok, {len(nodata)} returned no data, {len(failed)} FAILED')
    if nodata:
        print('\nreturned no data (legitimate on a small or partial run):')
        for name, _o, detail, _s in nodata:
            print(f'  {name}')
    if failed:
        print('\nFAILED -- these are the ones that matter:')
        for name, _o, detail, _s in failed:
            print(f'  {name}\n      {detail}')
        print('\nRe-run with --traceback for the full stack of each.')
    sys.exit(1 if failed else 0)


if __name__ == '__main__':
    main()
