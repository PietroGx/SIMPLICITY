#!/usr/bin/env python3
'''Trees, clustering, nextstrain and archiving, against the v2.4.71 layout.

    python tests/test_trees_clustering_archive.py

WHY

These four paths were untouched by the data flow refactor and untested after
it. All of them derive a destination path from a seeded output directory, and
that directory moved a level deeper (04_Output/<group>/<simulation>/<seed>)
with every component renamed. get_experiment_tree_simulation_dir,
get_clustering_table_filepath and get_nextstrain_dataset_paths all build their
targets from the SSOD helpers, and get_experiment_foldername_from_SSOD counted
back to path_parts[-4] until v2.4.71 -- which returns "04_Output" once a group
sits in between.

Self-contained: builds a real experiment with the serial runner and drives the
real functions over it. Nothing here needs the cluster.

archive_experiment DELETES the experiment folder it archives, so it runs last
and only ever against the scratch tree this test made.
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
import simplicity.output_manager as om
import simplicity.settings_manager as sm

failures = []
skipped = []


def check(label, got, want):
    ok = got == want
    print(f"  [{'ok ' if ok else 'FAIL'}] {label}: {got!r}")
    if not ok:
        failures.append(f'{label}: got {got!r}, want {want!r}')
    return ok


def check_runs(label, call, optional_deps=()):
    """Run it; a missing optional dependency is reported, not failed."""
    try:
        value = call()
    except ImportError as exc:
        print(f"  [skip] {label}: {exc}")
        skipped.append(f'{label}: {exc}')
        return None
    except Exception as exc:
        if any(d in str(exc) for d in optional_deps):
            print(f"  [skip] {label}: {type(exc).__name__}: {exc}")
            skipped.append(f'{label}: {exc}')
            return None
        print(f"  [FAIL] {label}: {type(exc).__name__}: {exc}")
        traceback.print_exc()
        failures.append(f'{label}: {type(exc).__name__}: {exc}')
        return None
    print(f"  [ok  ] {label}: {str(value)[:70]}")
    return value


def build(root, name):
    data_dir = os.path.join(root, 'Data')
    dm.set_data_dir(data_dir)
    sm._data_dir = data_dir
    import simplicity.runme as runme
    import simplicity.runners.serial as serial
    runme.run_experiment(
        name,
        lambda: ({'R': [1.2]},
                 {'population_size': 60, 'infected_individuals_at_start': 6,
                  'final_time': 25, 'long_shedders_ratio': 0.2,
                  'sequence_long_shedders': True}, 1),
        simplicity_runner=serial, archive_experiment=False)


def main():
    root = tempfile.mkdtemp(prefix='trees_probe_')
    name = 'trees_probe'
    try:
        print(f'scratch Data/ : {root}\nrunning a real experiment ...\n')
        build(root, name)

        sod = dm.get_simulation_output_dirs(name)[0]
        ssod = dm.get_seeded_simulation_output_dirs(sod)[0]
        print(f'\nagainst {os.path.relpath(ssod, root)}')

        print('\ntree destinations  (built from the SSOD helpers)')
        tree_dir = check_runs(
            'dm.get_experiment_tree_simulation_dir',
            lambda: dm.get_experiment_tree_simulation_dir(name, ssod))
        check_runs('dm.get_experiment_tree_simulation_files_dir',
                   lambda: dm.get_experiment_tree_simulation_files_dir(name, ssod))
        check_runs('dm.get_experiment_tree_simulation_plots_dir',
                   lambda: dm.get_experiment_tree_simulation_plots_dir(name, ssod))
        if tree_dir:
            # the tree directory is named for the SIMULATION, not the group and
            # not "04_Output" -- the path_parts[-4] bug landed exactly here
            check('tree dir is named for the simulation',
                  os.path.basename(tree_dir).startswith('sim_'), True)
            check('and it sits under 06_Trees',
                  dm.get_experiment_tree_dir(name) in tree_dir, True)

        # (experiment, ssod, tree_type, tree_subtype, file_type)
        for file_type in ('json', 'newick', 'img'):
            check_runs(f'om.get_tree_file_filepath ({file_type})',
                       lambda f=file_type: om.get_tree_file_filepath(
                           name, ssod, 'infection', 'full', f))
        check_runs('om.get_tree_filename',
                   lambda: om.get_tree_filename(name, ssod, 'infection',
                                                'full', 'json'))

        print('\nclustering')
        path = check_runs('om.get_clustering_table_filepath',
                          lambda: om.get_clustering_table_filepath(name, ssod))
        if path:
            import pandas as pd
            frame = pd.DataFrame({'lineage': ['wt'], 'cluster': [0]})
            check_runs('om.write_clustering_table',
                       lambda: om.write_clustering_table(name, ssod, frame)
                       or 'written')
            back = check_runs('om.read_clustering_table',
                              lambda: om.read_clustering_table(name, ssod))
            if back is not None:
                check('the clustering table round-trips', len(back), 1)
            check('the filename carries simulation and seed',
                  'seed_' in os.path.basename(path)
                  and 'sim_' in os.path.basename(path), True)

        print('\nnextstrain')
        check_runs('om.get_nextstrain_dataset_paths',
                   lambda: om.get_nextstrain_dataset_paths(name, ssod))

        print('\ntree building  (needs baltic 0.3.0 / ete3; skipped if absent)')
        check_runs('tree.build_tree.draw_tree (infection)',
                   lambda: _draw(name, ssod, 'infection'),
                   optional_deps=('baltic', 'ete3', 'ete4', 'Phylo'))

        print('\narchive_experiment  (tars the tree and DELETES the original)')
        experiment_dir = dm.get_experiment_dir(name)
        check_runs('om.archive_experiment',
                   lambda: om.archive_experiment(name) or 'archived')
        # it writes into Data/Archive/, not beside the experiment
        archive = os.path.join(dm.get_data_dir(), 'Archive', f'{name}.tar.gz')
        check('an archive was written', os.path.isfile(archive), True)
        check('and the original folder is gone',
              os.path.isdir(experiment_dir), False)
        if os.path.isfile(archive):
            import tarfile
            with tarfile.open(archive) as handle:
                members = handle.getnames()
            check('the archive contains the output tree',
                  any('/04_Output/main/' in m for m in members), True)
            # the two files that make an archived experiment interpretable at
            # all: what was run, and in what order. (state/ is Slurm-only --
            # the serial runner records no per-repeat state, so there is none
            # here to archive.)
            check('and the experiment record',
                  any(m.endswith('01_Settings/settings.json') for m in members),
                  True)
            check('and the repeat list',
                  any(m.endswith('03_Repeats/main/repeats.json') for m in members),
                  True)
    finally:
        shutil.rmtree(root, ignore_errors=True)

    print('\n' + '=' * 70)
    if skipped:
        print(f'{len(skipped)} skipped (optional dependency):')
        for line in skipped:
            print(f'  {line}')
    if failures:
        print(f'\n{len(failures)} FAILURE(S):')
        for line in failures:
            print(f'  {line}')
    else:
        print('\ntrees, clustering, nextstrain and archiving read the new layout')
    sys.exit(1 if failures else 0)


def _draw(name, ssod, tree_type):
    import simplicity.tree.build_tree as bt
    return bt.draw_tree(ssod, name, tree_type) or 'drawn'


if __name__ == '__main__':
    main()
