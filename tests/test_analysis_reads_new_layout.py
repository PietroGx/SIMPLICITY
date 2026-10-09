#!/usr/bin/env python3
'''The analysis path, against real output in the v2.4.71 layout.

    python tests/test_analysis_reads_new_layout.py

WHY THIS EXISTS

The data flow refactor is verified for setup, dispatch and the resolver, and a
real Slurm array runs correctly. What was never checked is the other half: the
40+ `get_parameter_value_from_simulation_output_dir` call sites, the
output_manager readers, the path helpers and the intra-host clock, all reading a
`04_Output/<group>/sim_NNN__<label>/seed_NNNN/` tree. Their signatures did not
change and they compile, which is not the same thing -- the tree they read moved
a level deeper and every directory was renamed.

So this runs a REAL experiment through runme.run_experiment with the serial
runner -- actual simulations, actual output files -- and then drives the actual
analysis functions over it. No stubs on the read side: the point is that these
exact functions, which the figures and the calibration fits call, work.

Small and short on purpose: this checks that the readers find and parse the
tree, not that any fit is well determined. A population of 60 over 20 days
cannot support a calibration regression and is not meant to -- and with few
diagnoses, some repeats have no sequencing output at all, which is correct
(see the note on sequencing below).

WHAT EACH GROUP OF CHECKS IS FOR

  path helpers     the SSOD helpers do path arithmetic, and the tree gained the
                   group level. get_experiment_foldername_from_SSOD counted back
                   to path_parts[-4] until v2.4.71, which returns "04_Output"
                   once a group sits between the output root and the simulation.

  readers          every output_manager read_* the figures and fits go through.

  resolver         the keystone, now against output a simulation actually wrote,
                   rather than directories a test created.

  intra-host clock evolutionary_rate.extract_ih_regression_data, which CLAUDE.md
                   names as the function diagnostics must reuse rather than
                   re-derive. If it cannot read the new tree, every calibration
                   is blocked.

  OSR table        write_OSR_vs_parameter_csv / read_OSR_vs_parameter_csv -- the
                   table both cal stages fit on. It walks the whole experiment,
                   so it exercises the listing and the resolver together.
'''
import os
import shutil
import sys
import tempfile
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import matplotlib
matplotlib.use('Agg')

import simplicity.dir_manager as dm
import simplicity.jobs as jobs
import simplicity.output_manager as om
import simplicity.settings_manager as sm
import simplicity.tuning.evolutionary_rate as er

failures = []


def check(label, got, want):
    ok = got == want
    print(f"  [{'ok ' if ok else 'FAIL'}] {label}: {got!r}")
    if not ok:
        failures.append(f'{label}: got {got!r}, want {want!r}')
    return ok


def check_runs(label, call):
    """For readers whose VALUE depends on the simulation: assert only that the
    real function reads the real tree without raising."""
    try:
        value = call()
    except Exception as exc:
        print(f"  [FAIL] {label}: {type(exc).__name__}: {exc}")
        traceback.print_exc()
        failures.append(f'{label}: {type(exc).__name__}: {exc}')
        return None
    summary = (f'{len(value)} rows' if hasattr(value, '__len__')
               else repr(value))
    print(f"  [ok  ] {label}: {summary}")
    return value


def build(root, name):
    """A real experiment: two simulations, two seeds, actually simulated."""
    data_dir = os.path.join(root, 'Data')
    dm.set_data_dir(data_dir)
    sm._data_dir = data_dir

    import simplicity.runme as runme
    import simplicity.runners.serial as serial
    runme.run_experiment(
        name,
        lambda: ({'R': [1.05, 1.3]},
                 {'population_size': 60,
                  'infected_individuals_at_start': 6,
                  'final_time': 20,
                  'long_shedders_ratio': 0.2,
                  'sequence_long_shedders': True}, 2),
        simplicity_runner=serial, archive_experiment=False)


def main():
    root = tempfile.mkdtemp(prefix='analysis_layout_')
    name = 'analysis_probe'
    try:
        print(f'scratch Data/ : {root}\nrunning a real experiment ...\n')
        build(root, name)

        sods = dm.get_simulation_output_dirs(name)
        ssods = [s for sod in sods
                 for s in dm.get_seeded_simulation_output_dirs(sod)]

        print('\nthe tree a simulation actually wrote')
        check('two simulations', len(sods), 2)
        check('four repeats', len(ssods), 4)
        check('under the group', [os.path.basename(
            os.path.dirname(os.path.dirname(s))) for s in ssods],
            ['main'] * 4)

        print('\npath helpers, against the group level they gained')
        for ssod in ssods[:1]:
            check('experiment name survives the extra level',
                  dm.get_experiment_foldername_from_SSOD(ssod), name)
            check('group is recoverable', dm.get_group_from_SSOD(ssod), 'main')
            check('simulation folder is the simulation, not the group',
                  dm.get_simulation_output_foldername_from_SSOD(ssod)
                  .startswith('sim_'), True)
            check('seed parses', dm.get_seed_from_SSOD(ssod), '0000')
        check('get_ssod finds a seed by number',
              os.path.basename(dm.get_ssod(sods[0], 1)), 'seed_0001')

        print('\nthe resolver, on output a simulation actually wrote')
        for sod in sods:
            value = sm.get_parameter_value_from_simulation_output_dir(sod, 'R')
            print(f'    {os.path.basename(sod)[:38]:<40} R = {value}')
        check('both swept values resolve exactly',
              {sm.get_parameter_value_from_simulation_output_dir(sod, 'R')
               for sod in sods}, {1.05, 1.3})

        print('\noutput_manager readers')
        ssod = ssods[0]
        check_runs('read_final_time', lambda: om.read_final_time(ssod))
        check_runs('read_individuals_data',
                   lambda: om.read_individuals_data(ssod))
        check_runs('read_phylogenetic_data',
                   lambda: om.read_phylogenetic_data(ssod))
        check_runs('read_lineage_frequency',
                   lambda: om.read_lineage_frequency(ssod))
        check_runs('read_simulation_trajectory',
                   lambda: om.read_simulation_trajectory(ssod))
        check_runs('read_sequencing_data_regression',
                   lambda: om.read_sequencing_data_regression(ssod))

        print('\nthe intra-host clock  (the function calibration must reuse)')
        check_runs('extract_ih_regression_data',
                   lambda: er.extract_ih_regression_data(ssod))
        check_runs('get_IH_lineages_data_experiment',
                   lambda: om.get_IH_lineages_data_experiment(name))

        print('\nthe OSR table both cal stages fit on')
        check_runs('write_OSR_vs_parameter_csv',
                   lambda: om.write_OSR_vs_parameter_csv(name, 'R'))
        # same min_seq_number / min_sim_lenght as the write: they are part of
        # the csv's filename, so a mismatch reads a file that was never written
        table = check_runs('read_OSR_vs_parameter_csv',
                           lambda: om.read_OSR_vs_parameter_csv(name, 'R', 0, 0))
        if table is not None and hasattr(table, 'columns'):
            check('the table carries the swept parameter',
                  'R' in table.columns, True)

        # The default read drops outliers, and detect_sod_outliers flags a lone
        # row (the statistic is degenerate on one observation), so with one
        # contributing repeat per simulation the default read is empty while
        # the csv is not. Assert the chain that is a layout question -- rows
        # written, rows read back, values right -- not the outlier policy.
        raw = check_runs('the csv itself has rows',
                         lambda: om.read_OSR_vs_parameter_csv(
                             name, 'R', 0, 0, include_outliers=True))
        if raw is not None and len(raw):
            check('one row per repeat that had sequencing output',
                  len(raw), sum(1 for one in ssods if os.path.isfile(
                      os.path.join(one, 'sequencing_data_regression.csv'))))
            check('and the swept values came through',
                  sorted(set(raw['R'])), [1.05, 1.3])

        # Row COUNT is not a layout question and must not be asserted here.
        # write_OSR_vs_parameter_csv wraps its per-repeat body in a bare
        # `except Exception: continue` (output_manager.py:627), and
        # er.tempest_regression raises on the single sequencing row a run this
        # small produces. An empty table therefore looks identical whether the
        # cause is thin data or a tree it cannot read -- so ask the layout
        # question directly, by doing per repeat exactly what that loop does.
        print('\n    ...so check the walk reaches every repeat, directly')
        always = ('read_final_time', 'read_individuals_data',
                  'read_phylogenetic_data', 'read_lineage_frequency',
                  'read_simulation_trajectory')
        reached, sequenced = [], []
        for one in ssods:
            try:
                for reader in always:
                    getattr(om, reader)(one)
                sm.get_parameter_value_from_simulation_output_dir(
                    os.path.dirname(one), 'R')
                reached.append(True)
            except Exception as exc:
                print(f'      {os.path.basename(one)}: {type(exc).__name__}: {exc}')
                reached.append(False)
            # A repeat with no DIAGNOSES has no sequencing file, and that is
            # correct, not a fault. There is no sequencing_rate any more: the
            # inline draw was replaced by recording t_diagnosis on every
            # diagnosed individual and writing the complete record, with any
            # subsample taken post-hoc (simplicity/sequencing.py). So the file
            # exists exactly when someone was diagnosed -- which is checked
            # below rather than assumed, because "no file" would otherwise look
            # the same as a tree the writer could not reach.
            individuals = om.read_individuals_data(one)
            diagnosed = (int(individuals['t_diagnosis'].notna().sum())
                         if 't_diagnosis' in individuals.columns else 0)
            sequenced.append((diagnosed, os.path.isfile(
                os.path.join(one, 'sequencing_data_regression.csv'))))

        check('every repeat is reachable and resolves', reached, [True] * 4)

        print('\n    the sequencing file tracks diagnoses, not a sampling rate')
        for one, (diagnosed, has_file) in zip(ssods, sequenced):
            print(f'      {os.path.basename(one)}: diagnosed={diagnosed} '
                  f'file={has_file}')
        check('a repeat has a sequencing file iff someone was diagnosed',
              [bool(d) == f for d, f in sequenced], [True] * 4)
        for one, (_d, has_file) in zip(ssods, sequenced):
            if has_file:
                check_runs(f'{os.path.basename(one)}: sequencing reads',
                           lambda o=one: om.read_sequencing_data_regression(o))
    finally:
        shutil.rmtree(root, ignore_errors=True)

    print('\n' + '=' * 70)
    if failures:
        print(f'{len(failures)} FAILURE(S) -- the analysis path does not read '
              f'the v2.4.71 layout:')
        for line in failures:
            print(f'  {line}')
    else:
        print('the analysis path reads the new layout')
    sys.exit(1 if failures else 0)


if __name__ == '__main__':
    main()
