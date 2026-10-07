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
Created on Tue Aug 27 19:38:43 2024

@author: pietro
"""

import os
import json
import simplicity.dir_manager as dm
import pandas as pd
import itertools
import copy

_data_dir = dm.get_data_dir()

def get_standard_parameters_values_file_path():
    standard_parameters_values_file_path = os.path.join(dm.get_reference_parameters_dir(), "standard_values.json")
    return standard_parameters_values_file_path

def get_parameter_specs_file_path():
    parameter_specs_file_path = os.path.join(dm.get_reference_parameters_dir(), "parameter_specs.json")
    return parameter_specs_file_path

_STANDARD_VALUES_DEFAULTS = {
        "population_size": 1000,
        'long_shedders_ratio': 0,
        "infected_individuals_at_start": 10,
        "tau_1": 2.86,
        "tau_2": 3.91,
        "tau_3": 7.5,
        "tau_3_long": 133.5,
        "tau_4": 8,
        "R": 1.1,
        "R_long": 1,
        "diagnosis_rate_standard": 0.1, # in percentage, will be converted to kds in model
        "diagnosis_rate_long"    : 0.1, # in percentage, will be converted to kdl in model
        "IH_virus_emergence_rate": 0,      # k_v in theoretical model equations
        "nucleotide_substitution_rate":  0.00008759,  # e in theoretical model equations
        "nucleotide_substitution_rate_long":  0.00008759,  # absolute NSR for long shedders
        "final_time": 365,
        "max_runtime": 86000,
        "phenotype_model": 'immune_waning',  # or 'linear'
        "consensus": 'argmax',  # or 'distribution' (immune_waning only)
        "sequence_long_shedders": False,
        "susceptibility_long": 1.0,
        "write_fasta": False,
        "seed": None
    }


def write_standard_parameters_values():
    filename= get_standard_parameters_values_file_path()
    standard_values = copy.deepcopy(_STANDARD_VALUES_DEFAULTS)
    with open(filename, "w") as file:
        json.dump(standard_values, file, indent=4)
    print(f"Standard values written to {filename}")

def write_parameter_specs():
    filename= get_parameter_specs_file_path()
    parameter_specs = {
        "population_size":               {"type": "int", "min": 0, "max": 10000},
        'long_shedders_ratio':           {"type": "float", "min": 0, "max": 1},
        "tau_1":                         {"type": "float", "min": 0, "max": 10},
        "tau_2":                         {"type": "float", "min": 0, "max": 100},
        "tau_3":                         {"type": "float", "min": 0, "max": 300},
        "tau_3_long":                    {"type": "float", "min": 0, "max": 400},
        "tau_4":                         {"type": "float", "min": 0, "max": 30},
        "infected_individuals_at_start": {"type": "int", "min": 0},
        "R":                             {"type": "float", "min": 0, "max": 20},
        "R_long":                        {"type": "float", "min": 0, "max": 20},
        "nucleotide_substitution_rate_long":  {"type": "float", "min": 0, "max": 1},
        "diagnosis_rate_standard":       {"type": "float", "min": 0, "max": 1},
        "diagnosis_rate_long":           {"type": "float", "min": 0, "max": 1},
        "IH_virus_emergence_rate":       {"type": "float", "min": 0},
        "nucleotide_substitution_rate":  {"type": "float", "min": 0, "max": 1},
        "final_time":                    {"type": "int", "min": 0},    
        "max_runtime":                   {"type": "int", "min": 0},
        "phenotype_model":               {"type": "str"},
        "consensus":                     {"type": "str"},
        "sequence_long_shedders":        {"type": "bool"},
        "susceptibility_long":           {"type": "float", "min": 0},
        "write_fasta":                   {"type": "bool"}
        }

    with open(filename, "w") as file:
        json.dump(parameter_specs, file, indent=4)
    print(f"Parameter specifications written to {filename}")

def read_standard_parameters_values():
    filename= get_standard_parameters_values_file_path()
    try:
        with open(filename, "r") as file:
            from_file = json.load(file)
    except FileNotFoundError:
        print(f"Error: {filename} not found. Writing default standard values.")
        write_standard_parameters_values()
        with open(filename, "r") as file:
            from_file = json.load(file)

    # Code defaults under the file, so a parameter added after this file was
    # written still resolves instead of raising KeyError downstream.
    return {**_STANDARD_VALUES_DEFAULTS, **from_file}

def read_parameter_specs():
    filename= get_parameter_specs_file_path()
    try:
        with open(filename, "r") as file:
            return json.load(file)
    except FileNotFoundError:
        print(f"Error: {filename} not found. Writing default parameter specifications.")
        write_parameter_specs()
        return read_parameter_specs()
    
def write_user_set_parameters_file(user_set_parameters, filename):
    file_path = os.path.join(dm.get_reference_parameters_dir(),filename)
    with open(file_path, "w") as file:
        json.dump(user_set_parameters, file, indent=4)
    print(f"user_set_parameters saved to {file_path}")

def read_user_set_parameters_file(filename):
    file_path = os.path.join(dm.get_reference_parameters_dir(),filename)
    try:
        with open(file_path, "r") as file:
            return json.load(file)
    except FileNotFoundError:
        print(f"Error: {filename} not found. Writing default standard values.")
        write_standard_parameters_values()
        return read_standard_parameters_values()
    
# A _scenario_groups entry may name itself with this key; unnamed groups get a
# positional name. A plain experiment is one group, DEFAULT_GROUP.
GROUP_NAME_KEY = 'name'
DEFAULT_GROUP = 'main'


def get_experiment_settings_file_path(experiment_name):
    """The experiment record: n_seeds, its groups, and its ordered simulations.
    One file -- it used to be two, written together and read apart, with n_seeds
    read by two separate functions that did the same thing."""
    return os.path.join(dm.get_experiment_settings_dir(experiment_name),
                        'settings.json')

def check_parameters_names(parameters_dic):
    STANDARD_VALUES = read_standard_parameters_values()
    for key in parameters_dic.keys():
        if key not in STANDARD_VALUES.keys():
            raise ValueError(f'Parameter {key} is not a valid parameter')

def read_experiment_settings_file(experiment_name):
    """{'n_seeds', 'groups', 'simulations'} -- the whole record."""
    with open(get_experiment_settings_file_path(experiment_name)) as json_file:
        return json.load(json_file)


def read_simulations(experiment_name):
    """The ordered simulation records: {'id', 'group', 'parameters'}."""
    return read_experiment_settings_file(experiment_name)['simulations']


def read_groups(experiment_name):
    """The ordered group records: {'name', 'n_seeds'}."""
    return read_experiment_settings_file(experiment_name)['groups']


def read_experiment_settings(experiment_name):
    """Just the parameter dicts, in order -- the pre-groups contract."""
    return [record['parameters'] for record in read_simulations(experiment_name)]


def read_n_seeds_file(experiment_name):
    """Returns a dict, so the existing ['n_seeds'] callers are unchanged."""
    return {'n_seeds': read_experiment_settings_file(experiment_name)['n_seeds']}


def _as_simulations(tagged):
    """Attach the identity. The id is the position in a deterministic order
    (itertools.product over declared keys, groups in declared order), so it is
    reproducible and cannot collide. It stays OUTSIDE 'parameters':
    check_parameters_names rejects any key that is not a real parameter."""
    return [{'id': index, 'group': group, 'parameters': parameters}
            for index, (group, parameters) in enumerate(tagged)]


def generate_experiment_settings(varying_params: dict, fixed_params: dict = None):
    """
    Generates a list of parameter combinations from varying and fixed parameters.

    Args:
        varying (dict): Parameters for which all combinations should be generated.
                        An optional reserved key '_scenario_groups' may hold a list
                        of dicts, one per "group" (e.g. one per experiment scenario).
                        Within each group dict, list/tuple-valued entries are that
                        group's own local varying axis (expanded independently of
                        every other group); scalar entries are fixed overrides for
                        that group only. This lets correlated parameter sets (e.g.
                        different NSR sweep ranges per scenario) be submitted as a
                        single combined experiment instead of one experiment per
                        group. A group may carry a 'name' key to name itself;
                        unnamed groups are named by position.
        fixed (dict): Parameters that should have the same value across all combinations.

    Returns:
        List[dict]: ordered simulation records {'id', 'group', 'parameters'}.
                    The grouping used to be flattened away here, which is why
                    anything needing it downstream either re-invoked the config
                    or re-invented a label at the plot call site.
    """
    fixed_params = fixed_params or {}
    varying_params = dict(varying_params or {})
    scenario_groups = varying_params.pop('_scenario_groups', None)

    if scenario_groups is None:
        keys, values = zip(*varying_params.items()) if varying_params else ([], [])
        combinations = list(itertools.product(*values)) if values else [()]

        experiment_settings = []
        for combo in combinations:
            setting = dict(zip(keys, combo))
            setting.update(copy.deepcopy(fixed_params))  # Avoid mutation
            experiment_settings.append((DEFAULT_GROUP, setting))

        return _as_simulations(experiment_settings)

    # Grouped path: each group expands its own list-valued keys independently,
    # so different groups can vary different parameters over different ranges.
    experiment_settings = []
    for group_index, group in enumerate(scenario_groups):
        group_name = group.get(GROUP_NAME_KEY) or f'group_{group_index:02d}'
        # the name is identity, not a parameter: keep it out of the spec
        group_spec = {k: v for k, v in group.items() if k != GROUP_NAME_KEY}
        group_varying = {k: v for k, v in group_spec.items() if isinstance(v, (list, tuple))}
        group_fixed = {k: v for k, v in group_spec.items() if not isinstance(v, (list, tuple))}

        keys, values = zip(*group_varying.items()) if group_varying else ([], [])
        combinations = list(itertools.product(*values)) if values else [()]

        for combo in combinations:
            setting = dict(zip(keys, combo))
            setting.update(copy.deepcopy(fixed_params))
            setting.update(copy.deepcopy(group_fixed))
            experiment_settings.append((group_name, setting))

    return _as_simulations(experiment_settings)

def write_experiment_settings(experiment_name: str, experiment_settings: list, n_seeds: int):
    """
    Writes experiment settings (a list of parameter dictionaries) to a JSON file.

    Args:
        experiment_name (str): Name of the experiment (used for output folder).
        experiment_settings (list): List of parameter dictionaries.
        n_seeds (int): Number of random seeds to be stored separately.
    """
    # check parameter names validity
    for record in experiment_settings:
        check_parameters_names(record['parameters'])

    # One group record per distinct group, in first-seen order. n_seeds is per
    # group: the pipeline drivers already carry --cal-seeds and --exp-seeds as
    # separate knobs, which only worked while each stage was its own experiment.
    groups = [{'name': name, 'n_seeds': n_seeds}
              for name in dict.fromkeys(r['group'] for r in experiment_settings)]

    experiment_settings_file_path = get_experiment_settings_file_path(experiment_name)
    with open(experiment_settings_file_path, 'w') as settings_file:
        json.dump({'n_seeds': n_seeds,
                   'groups': groups,
                   'simulations': experiment_settings}, settings_file, indent=4)

    print(f"Experiment settings file written to {experiment_settings_file_path}")

def write_simulation_parameters(file_path, settings):
    """One simulation's complete parameters.

    Takes a dict. This was 23 positional arguments threaded from a dict at the
    only call site, where a transposed pair of same-typed values was silent and
    every new parameter meant editing a signature, a call and a dict in step.
    """
    with open(file_path, "w") as json_file:
        json.dump({**settings, "t_0": 0}, json_file, indent=4)


def generate_filename_from_params(params: dict):
    
    abbreviations = {
    "population_size": "N",
    "tau_3": "tau3",
    "infected_individuals_at_start": "init",
    "R": "R",
    "R_long": "Rl",
    "nucleotide_substitution_rate_long": "NSR_long",
    "diagnosis_rate_standard": "kds",
    "diagnosis_rate_long": "kdl",
    "IH_virus_emergence_rate": "kv",
    "nucleotide_substitution_rate": "NSR",
    "final_time": "T",
    "phenotype_model": "pheno",
    "consensus": "cdist"
    # excluded: max_runtime, seed, F
}
    exclude  = {"max_runtime", "seed"}
    parts = []
    for key, value in params.items():
        if key in exclude:
            continue
        abbrev = abbreviations.get(key, key)
        if isinstance(value, float):
            # 6 significant figures, not 2. Nothing resolves through this string
            # any more -- simulation_stem's id prefix carries identity and
            # parameters.json carries the values -- but people still read it off
            # `ls`, and at 2 figures a folder named R_1p1 held R=1.06.
            value_str = f"{value:.6g}".replace('.', 'p')
        elif isinstance(value, int):
            value_str = str(value)
        elif isinstance(value, str):
            value_str = value.replace(' ', '')
        else:
            value_str = str(value)
        parts.append(f"{abbrev}_{value_str}")
    
    file_name = "_".join(parts) + ".json"
    return file_name


# One simulation's parameters, written beside the output they describe. This is
# what get_parameter_value_from_simulation_output_dir reads, so nothing has to
# rebuild a filename to find them.
PARAMETERS_FILE = 'parameters.json'

# Keep a stem well inside the 255-byte filename limit; the id makes it unique,
# so the label is what gets cut.
MAX_STEM = 180


def simulation_stem(record, standard_values=None):
    """Directory and file stem for one simulation: its id, then a label.

    The id is the identity -- assigned, dense, reproducible, and impossible to
    collide. The label is for humans. Two parameter sets whose labels coincide
    used to share one filename and silently overwrite each other; now they
    cannot, because the stems differ in the id.
    """
    standard_values = standard_values or read_standard_parameters_values()
    modified = {key: value for key, value in record['parameters'].items()
                if key in standard_values and value != standard_values[key]}
    label = generate_filename_from_params(modified)[:-len('.json')] if modified \
        else 'standard_values'
    prefix = f"sim_{record['id']:03d}__"
    return prefix + label[:MAX_STEM - len(prefix)]


def read_settings_and_write_simulation_parameters(experiment_name):
    """
    Reads an experiment settings file in JSON format and generates individual simulation parameter files
    based on the combinations of parameters in the settings. The function will create a separate JSON file
    for each parameter combination within a directory named after the experiment.

    Parameters:
    -----------
    experiment_name : str
        The name of the experiment. This is used to locate the settings file and to create the corresponding
        simulation parameters directory.

    """
    
    STANDARD_VALUES = read_standard_parameters_values()
    parameters_dir = dm.get_simulation_parameters_dir(experiment_name)
    experiment_output_dir = dm.get_experiment_output_dir(experiment_name)

    for record in read_simulations(experiment_name):
        stem = simulation_stem(record, STANDARD_VALUES)
        settings = {**STANDARD_VALUES, **record['parameters']}

        write_simulation_parameters(
            os.path.join(parameters_dir, f'{stem}.json'), settings)

        # The reader-facing copy, written here at setup rather than by the task:
        # every repeat of this simulation would otherwise race to write it.
        # Under the group, so a group's output subtree is self-contained.
        simulation_output_dir = os.path.join(experiment_output_dir,
                                             record['group'], stem)
        os.makedirs(simulation_output_dir, exist_ok=True)
        write_simulation_parameters(
            os.path.join(simulation_output_dir, PARAMETERS_FILE), settings)

    print(f"Simulation parameters written to directory: {parameters_dir}")

def read_simulation_parameters(experiment_name, stem):
    """One simulation's complete parameters, by the stem repeats.json records.

    One small file, which is why 02_Simulations is kept: a task would otherwise
    parse settings.json, holding every group's parameters, to find its own.
    """
    with open(os.path.join(dm.get_simulation_parameters_dir(experiment_name),
                           f'{stem}.json')) as handle:
        return json.load(handle)


def get_simulation_parameters(simulation_output_dir):
    """The parameters that produced this output, read from beside it."""
    with open(os.path.join(simulation_output_dir, PARAMETERS_FILE)) as file:
        return json.load(file)


def get_parameter_value_from_simulation_output_dir(simulation_output_dir, parameter):
    """Read one parameter of the simulation that produced this output.

    Same signature as ever -- 43 call sites depend on it. What changed is that
    it no longer rebuilds 02_Simulations/<folder name>.json from
    parts[-3] and parts[-1] of the path it was handed. That string coupling made
    the folder name a lookup key, and the folder name is a lossy encoding of the
    parameters, so two swept values that rounded together resolved to one file.
    """
    return get_simulation_parameters(simulation_output_dir)[parameter]


def read_OSR_NSR_regressor_parameters():
    file_path = os.path.join(dm.get_reference_parameters_dir(),
              'OSR_NSR_regressor_parameters_for_standard_parameter_values_exp.csv')
    df = pd.read_csv(file_path,index_col=0)
    best_fit_df = pd.to_numeric(df['Best Fit'], errors='coerce')
    return best_fit_df
    
