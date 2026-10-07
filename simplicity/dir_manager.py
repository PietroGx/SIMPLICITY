# This file is part of SIMPLICITY
# Copyright (C) 2025 Pietro Gerletti
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 19:36:24 2024

@author: pietro
"""

import os

# Default data directory 
DEFAULT_DATA_DIR = os.path.join(os.getcwd(), "Data")
os.makedirs(DEFAULT_DATA_DIR,exist_ok=True)

# Global variable to store the data directory path
_data_dir = DEFAULT_DATA_DIR

# set env variables
os.environ["SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_MULTIPROCESS"] = str(10)
os.environ["SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM"] = str(200)
os.environ["QT_QPA_PLATFORM"] = "offscreen"
# Per-task SLURM resource request defaults (2G/2 days). Memory is from real
# sacct data on the isolated long-NSR calibration grid -- see
# impact_long_shedders_config. Time was bumped from 1 day after that same
# grid showed some high-R_long tau groups (e.g. edge_case) still running
# past the old 1-day walltime and getting killed mid-simulation (see
# simplicity.runners.slurm.reconcile_terminated_tasks, which now detects and
# signals that case regardless, but avoiding the kill in the first place is
# still better when a day's headroom is cheap to give).
# CLI scripts that expose --slurm-mem/--slurm-time overwrite these with the
# user's chosen value before submitting; this is just the baseline.
os.environ["SIMPLICITY_SLURM_MEM"] = "4G"
os.environ["SIMPLICITY_SLURM_TIME"] = "2-00:00:00"

def set_data_dir(path):
    """Set the data directory path."""
    global _data_dir
    _data_dir = path

def get_data_dir():
    """Get the current data directory path."""
    return _data_dir

def get_reference_parameters_dir():
    reference_parameters_dir = os.path.join(get_data_dir(), '00_Reference_parameters')
    os.makedirs(reference_parameters_dir, exist_ok=True)
    return reference_parameters_dir

# The numbered subdirectories, named once. They were spelled as literals in
# eleven places, so a rename meant finding all eleven.
SETTINGS_DIRNAME    = "01_Settings"
SIMULATIONS_DIRNAME = "02_Simulations"
REPEATS_DIRNAME     = "03_Repeats"
OUTPUT_DIRNAME      = "04_Output"


def create_directories(experiment_name):
    """Create necessary subdirectories within the data directory."""
    
    experiment_dir = os.path.join(_data_dir,experiment_name)
    subdirectories = [SETTINGS_DIRNAME,
                      SIMULATIONS_DIRNAME,
                      REPEATS_DIRNAME,
                      OUTPUT_DIRNAME,
                      ]
    
    for subdir in subdirectories:
        path = os.path.join(experiment_dir, subdir)
        os.makedirs(path, exist_ok=True)

def get_experiment_dir(experiment_name):
    """Get the experiment_name directory path."""
    experiment_dir = os.path.join(_data_dir,experiment_name)
    if os.path.isdir(experiment_dir):
        return experiment_dir
    else:
        raise ValueError('No experiment with that name!')
        
def get_experiment_settings_dir(experiment_name):
    """Get the experiment_name simulation parameters directory path."""
    experiment_dir = os.path.join(_data_dir,experiment_name)
    if os.path.isdir(experiment_dir):
        return os.path.join(experiment_dir, SETTINGS_DIRNAME)
    else:
         raise ValueError('No experiment with that name!')

def get_simulation_parameters_dir(experiment_name):
    """Get the experiment_name simulation parameters directory path."""
    experiment_dir = os.path.join(_data_dir,experiment_name)
    if os.path.isdir(experiment_dir):
        return os.path.join(experiment_dir, SIMULATIONS_DIRNAME)
    else:
         raise ValueError('No experiment with that name!')
         
def get_repeats_dir(experiment_name, group=None):
    """The repeats directory: one ordered repeat list per group, plus that
    group's state and progress. Replaces 03_Seeded_simulation_parameters, which
    held a copy of the combination's parameters per seed plus five signal files
    that carried one bit each by existing."""
    base = os.path.join(get_experiment_dir(experiment_name), REPEATS_DIRNAME)
    path = base if group is None else os.path.join(base, group)
    os.makedirs(path, exist_ok=True)
    return path


def get_experiment_output_dir(experiment_name):
    """Get the experiment_name output directory path."""
    experiment_dir = os.path.join(_data_dir,experiment_name)
    if os.path.isdir(experiment_dir):
        return os.path.join(experiment_dir, OUTPUT_DIRNAME)
    else:
         raise ValueError('No experiment with that name!')
         
def get_experiment_plots_dir(experiment_name):
    experiment_dir = os.path.join(_data_dir,experiment_name)
    experiment_plots_dir = os.path.join(experiment_dir, '05_Plots')
    os.makedirs(experiment_plots_dir, exist_ok=True) 
    return experiment_plots_dir

def get_experiment_simulations_plots_dir(experiment_name):
    experiment_plots_dir = get_experiment_plots_dir(experiment_name)
    experiment_simulations_plots_dir = os.path.join(experiment_plots_dir, 'Simulations')
    os.makedirs(experiment_simulations_plots_dir, exist_ok=True) 
    return experiment_simulations_plots_dir

def get_experiment_tree_dir(experiment_name):
    experiment_dir = get_experiment_dir(experiment_name)
    experiment_tree_dir = os.path.join(experiment_dir, '06_Trees')
    os.makedirs(experiment_tree_dir, exist_ok=True) 
    return experiment_tree_dir

def get_nextstrain_dir(experiment_name):
    base = get_experiment_tree_dir(experiment_name)
    ns_dir = os.path.join(base, "nextstrain")
    os.makedirs(ns_dir, exist_ok=True)
    return ns_dir

def get_experiment_cluster_dir(experiment_name):
    exp_dir = get_experiment_dir(experiment_name)
    out = os.path.join(exp_dir, "07_Clustering")
    os.makedirs(out, exist_ok=True)
    return out

def get_experiment_fit_result_dir(experiment_name):
    fit_result_dir = os.path.join(get_experiment_dir(experiment_name), 'Fit_results')
    os.makedirs(fit_result_dir, exist_ok=True)
    return fit_result_dir

def get_slurm_logs_dir(experiment_name):
    slurm_logs_dir = os.path.join(get_experiment_dir(experiment_name), 'slurm','slurm_logs')
    os.makedirs(slurm_logs_dir, exist_ok=True)
    return slurm_logs_dir

def get_slurm_id_map_dir(experiment_name):
    slurm_id_map_dir  = os.path.join(get_experiment_dir(experiment_name), 'slurm','job_id_mapping')
    os.makedirs(slurm_id_map_dir , exist_ok=True)
    return slurm_id_map_dir

def get_simulation_output_dirs(experiment_name, group=None):
    """Every simulation's output directory, or only one group's.

    group=None returns all of them, which is what every existing caller gets.
    Passing a group is how a read says which part of a heterogeneous experiment
    it means: an experiment already holds several independent sweeps (one
    _scenario_groups entry each), and a listing that mixes them feeds a plot
    series sorting by parameter value across sweeps it should not compare.

    Ordered: groups as declared, then simulations by id (the sim_NNN prefix
    makes the string sort an id sort, up to 999 simulations).
    """
    experiment_output_dir = get_experiment_output_dir(experiment_name)
    if group is None:
        groups = get_groups(experiment_name)
    else:
        if group not in get_groups(experiment_name):
            raise ValueError(f'{experiment_name} has no group {group!r}; it '
                             f'has {get_groups(experiment_name)}')
        groups = [group]

    simulation_output_dirs = []
    for name in groups:
        group_dir = os.path.join(experiment_output_dir, name)
        if not os.path.isdir(group_dir):
            continue
        simulation_output_dirs.extend(sorted(
            os.path.join(group_dir, folder)
            for folder in os.listdir(group_dir)
            if os.path.isdir(os.path.join(group_dir, folder))))
    return simulation_output_dirs


def get_group_from_SSOD(seeded_simulation_output_dir):
    """The group of Data/<exp>/04_Output/<group>/<simulation>/<seed>."""
    return os.path.basename(os.path.dirname(
        os.path.dirname(seeded_simulation_output_dir)))


def get_groups(experiment_name):
    """The experiment's group names, in declared order."""
    import simplicity.settings_manager as sm
    return [g['name'] for g in sm.read_groups(experiment_name)]


def get_seeded_simulation_output_dirs(simulation_output_dir):
    """This simulation's repeats, sorted -- so seed_0000 is always first."""
    return sorted(
        os.path.join(simulation_output_dir, subfolder)
        for subfolder in os.listdir(simulation_output_dir)
        if os.path.isdir(os.path.join(simulation_output_dir, subfolder)))

# SSOD = seeded_simulation_output_dir
def get_simulation_output_foldername_from_SSOD(seeded_simulation_output_dir):
    ''' Get the simulation_output folder_name of Data/experiment/04_Output/simulation_output/seed_nr'''
    simulation_output_foldername = os.path.basename(os.path.dirname(seeded_simulation_output_dir))
    return simulation_output_foldername

def get_experiment_foldername_from_SSOD(seeded_simulation_output_dir):
    """The experiment folder name, found by locating the output directory
    component rather than counting back a fixed number of parts.

    It was path_parts[-4], which is only right while the tree is exactly
    <experiment>/04_Output/<simulation>/<seed>. Add or remove one level and it
    silently returns a different component -- "04_Output" as the experiment
    name, say -- with nothing to signal it.
    """
    path_parts = os.path.normpath(seeded_simulation_output_dir).split(os.sep)
    for index in range(len(path_parts) - 1, 0, -1):     # last match wins
        if path_parts[index] == OUTPUT_DIRNAME:
            return path_parts[index - 1]
    raise ValueError(f'Not a path under {OUTPUT_DIRNAME}/: '
                     f'{seeded_simulation_output_dir}')

def get_experiment_tree_simulation_dir(experiment_name,
                                       seeded_simulation_output_dir):
    experiment_tree_dir = get_experiment_tree_dir(experiment_name)
    simulation_output_foldername = get_simulation_output_foldername_from_SSOD(seeded_simulation_output_dir)
    experiment_tree_simulation_dir = os.path.join(experiment_tree_dir,simulation_output_foldername)
    os.makedirs(experiment_tree_simulation_dir,exist_ok=True)
    return experiment_tree_simulation_dir

def get_experiment_tree_simulation_files_dir(experiment_name,
                                       seeded_simulation_output_dir):
    experiment_tree_simulation_dir = get_experiment_tree_simulation_dir(experiment_name,
                                           seeded_simulation_output_dir)
   
    experiment_tree_simulation_files_dir = os.path.join(experiment_tree_simulation_dir,
                                                       'files')
    os.makedirs(experiment_tree_simulation_files_dir,exist_ok=True)
    return experiment_tree_simulation_files_dir

def get_experiment_tree_simulation_plots_dir(experiment_name,
                                       seeded_simulation_output_dir):
    experiment_tree_simulation_dir = get_experiment_tree_simulation_dir(experiment_name,
                                           seeded_simulation_output_dir)
    experiment_tree_simulation_plots_dir = os.path.join(experiment_tree_simulation_dir,
                                                       'plots')
    os.makedirs(experiment_tree_simulation_plots_dir,exist_ok=True)
    return experiment_tree_simulation_plots_dir

def get_ssod(sim_out_dir, seed_number):
    """
    Returns the seeded simulation output directory (SSOD) for the given seed number.
    """
    all_seed_dirs = get_seeded_simulation_output_dirs(sim_out_dir)
    target = f"seed_{seed_number:04d}"
    
    for path in all_seed_dirs:
        if os.path.basename(path) == target:
            return path

    raise ValueError(f"Seed folder '{target}' not found in {sim_out_dir}")
    
def get_seed_from_SSOD(seeded_simulation_output_dir):
    """
    Extract the zero-padded seed number (e.g. '0007')
    """
    import re
    m = re.search(r"seed_(\d+)$", seeded_simulation_output_dir)
    if not m:
        raise ValueError(f"Could not find 'seed_####' at the end of: {seeded_simulation_output_dir}")
    return m.group(1)

# Outputs a run may be told not to write, by base name without ".csv". Set
# SIMPLICITY_SKIP_OUTPUTS to a comma-separated list. Lives here because both
# population (which streams two of them) and output_manager (which saves them)
# need it, and importing one from the other closes a cycle.
#
# Why this exists: a calibration run reads only final_time.csv,
# sequencing_data_regression.csv, individuals_data.csv and
# phylogenetic_data.csv -- the first two for the OSR fit, the other two through
# evolutionary_rate.extract_ih_regression_data for the intra-host clock. It
# wrote three more files that nothing in the calibration or sanity path ever
# opened. Measured on the calibrated grid: lineage_frequency.csv alone reached
# 100.6 GB, one simulation of the top NSR_long sweep point weighing 1.1 GB --
# a grid point that exists to be rejected by the fit.
SKIP_OUTPUTS_ENV = "SIMPLICITY_SKIP_OUTPUTS"


def skipped_outputs():
    """Base names of the outputs this run must not write."""
    raw = os.environ.get(SKIP_OUTPUTS_ENV, "")
    return frozenset(name.strip() for name in raw.split(",") if name.strip())


def skipping(name):
    return name in skipped_outputs()
