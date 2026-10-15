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
@author: pietro

STANDARD_VALUES for SIMPLICITY simulation: 
    
    "population_size": 1000
    "infected_individuals_at_start": 100
    "R": 1.5
    "k_d": 0.0055
    "k_v": 0.0085
    "e": 0.0017 (evolutionary rate)
    "final_time": 365*3 
    "max_runtime": 300
    "phenotype_model": 'immune waning' or 'distance from wt'
    "seed": None
    "F": 1.25
    
If you want to change any, you can specify them in the parameters dictionary below. 
For each parameter, specify a list of values that you would like to use for the 
simulation. If you want to change more than one parameter at the time, consider 
that you need to enter the same number of values for each parameter, e.g. :
    par 1 = [value1, value2]
    par 2 = [value3, value4]
This will run a simulation with par 1 = value1 and par 2 = value 3, and a simulation
with par 1 = value2 and par 2 = value4. 

Each simulation will be repeated n_seeds time with a different random seed.

The set of all simulations is what we call an experiment.
"""
import simplicity.dir_manager as dm
import simplicity.jobs            as jobs
import simplicity.settings_manager as sm
import simplicity.output_manager as om
import simplicity.runners.serial 
from   simplicity.runners.unit_run import run_seeded_simulation
import types

experiment_name = 'Test_runme'
# experiment settings
def user_set_experiment_settings():
    varying_params = {
        'phenotype_model': ['distance from wt', 'immune waning']
    }
    # parameters to keep fixed (but different from standard_value) across combinations
    fixed_params = {}
    
    n_seeds         = 2
    return (varying_params,fixed_params,n_seeds)

def run_experiment_phase(experiment_name: str,
                         set_experiment_parameters: types.FunctionType,
                         simplicity_runner: types.ModuleType,
                         append: bool = True):
    """Add one PHASE's groups to an experiment and run only those groups.

    Stage 8's model: an experiment is the whole study, and a phase -- a
    calibration sweep, or the experiment-phase scenarios -- is a set of named
    groups inside it. The phases are sequential and data-dependent (the
    experiment phase's parameters come from the calibration's fit), so they
    cannot all be declared up front; each is appended when its inputs exist.

    Differs from run_experiment in exactly two ways, both of which have to
    hold together:

      the simulation ids CONTINUE from what is already there (next_simulation_id),
      because simulation_stem puts the id in the output directory name and a
      second phase restarting at 0 would overwrite the first one's output;

      the runner is given THIS phase's groups, so appending does not re-run
      the groups that are already finished.

    append=False writes the experiment from scratch, which is what the first
    phase of a run wants.
    """
    dm.create_directories(experiment_name)
    varying_params, fixed_params, n_seeds = set_experiment_parameters()

    start_id = sm.next_simulation_id(experiment_name) if append else 0
    experiment_settings = sm.generate_experiment_settings(
        varying_params, fixed_params, start_id=start_id)

    # the groups this phase is adding, in declared order and without
    # duplicates -- what the runner is scoped to below
    phase_groups = list(dict.fromkeys(record['group']
                                      for record in experiment_settings))

    sm.write_experiment_settings(experiment_name, experiment_settings, n_seeds,
                                 append=append)
    sm.read_settings_and_write_simulation_parameters(experiment_name)
    # per group, and it refuses to renumber a group that already carries state
    jobs.write_repeats(experiment_name)

    print('')
    print('-------------------------------------------------------------------')
    print(f'{experiment_name}: running phase groups {phase_groups}')
    simplicity_runner.run_seeded_simulations(experiment_name,
                                             run_seeded_simulation,
                                             groups=phase_groups)
    print('')
    print('-------------------------------------------------------------------')
    print(f'PHASE {phase_groups} OF {experiment_name} EXECUTED SUCCESSFULLY.')
    print('###################################################################')
    return phase_groups


def run_experiment(experiment_name: str, 
                   set_experiment_parameters: types.FunctionType,
                   simplicity_runner: types.ModuleType,
                   archive_experiment=False):
    
    # setup experiment files directories
    dm.create_directories(experiment_name)
    # set parameters 
    varying_params, fixed_params, n_seeds = set_experiment_parameters()
    experiment_settings = sm.generate_experiment_settings(varying_params, fixed_params)
    # Write experiment settings file
    sm.write_experiment_settings(experiment_name, experiment_settings, n_seeds)
    # write simulation parameters files
    sm.read_settings_and_write_simulation_parameters(experiment_name)
    # write the ordered repeat list: one entry per (simulation, seed), which is
    # what the runners address and what the Slurm array maps onto
    jobs.write_repeats(experiment_name)
    print('')
    print('-------------------------------------------------------------------')
    # let one of simplicity.runners run each seeded simulation
    simplicity_runner.run_seeded_simulations(experiment_name, 
                                             run_seeded_simulation)
    
    if archive_experiment: 
        om.archive_experiment(experiment_name)
    
    print('')
    print('-------------------------------------------------------------------')
    print(f'EXPERIMENT {experiment_name} EXECUTED SUCCESSFULLY.')
    print('###################################################################')

if __name__ == "__main__":
    
    run_experiment(experiment_name,       
                   user_set_experiment_settings,
                   simplicity_runner  = simplicity.runners.serial,
                   archive_experiment = False)