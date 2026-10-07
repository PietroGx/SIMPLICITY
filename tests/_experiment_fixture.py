#!/usr/bin/env python3
'''A minimal real experiment on disk, for tests of the Slurm control plane.

Built through the actual settings and jobs path -- no simulation runs -- so a
test exercises the same files the pipeline writes rather than a stand-in.

These tests used to monkeypatch sm.get_seeded_simulation_parameters_paths and
touch signal files by hand. That made them pass against a shape the pipeline no
longer had: a diagnostic that reimplements its stage stops describing that stage
the moment the stage changes, which is the failure CLAUDE.md is about.
'''
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import simplicity.dir_manager as dm
import simplicity.jobs as jobs
import simplicity.settings_manager as sm

# run_seeded_simulations reads this with a hard subscript
os.environ.setdefault('SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM', '4')


class Experiment:
    '''A real experiment in a scratch Data/, with repeats but no output.'''

    def __init__(self, n_simulations=4, n_seeds=1, name='probe', groups=None):
        self.root = tempfile.mkdtemp(prefix='slurm_test_')
        # both, deliberately: sm._data_dir is a copy taken at import time and
        # does not follow dm.set_data_dir
        dm.set_data_dir(self.root)
        sm._data_dir = self.root
        self.name = name
        dm.create_directories(name)

        if groups is None:
            varying = {'nucleotide_substitution_rate':
                       [1e-4 * (i + 1) for i in range(n_simulations)]}
        else:
            varying = {'_scenario_groups': groups}
        sm.write_experiment_settings(
            name,
            sm.generate_experiment_settings(varying,
                                            {'population_size': 20,
                                             'final_time': 2}),
            n_seeds)
        sm.read_settings_and_write_simulation_parameters(name)
        jobs.write_repeats(name)
        # the submission map's shape, which is what the control plane iterates
        self.pairs = [[group, record['index']]
                      for group, record in jobs.all_repeats(name)]

    def key(self, position):
        group, index = self.pairs[position]
        return group, index

    def set(self, position, state, **extra):
        group, index = self.key(position)
        return jobs.set_state(self.name, group, index, state, **extra)

    def state(self, position):
        group, index = self.key(position)
        return jobs.state_of(self.name, group, index)

    def states(self):
        return [self.state(i) for i in range(len(self.pairs))]

    def close(self):
        shutil.rmtree(self.root, ignore_errors=True)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
