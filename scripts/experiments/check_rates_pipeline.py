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
   
Each simulation will be repeated n_seeds time with a different random seed.

The set of all simulations is what we call an experiment.
"""

from experiment_script_runner import run_experiment_script
import simplicity.settings_manager as sm
import simplicity.plots_manager as pm
import argparse

experiment_name = 'check_rates_pipeline'

def user_set_experiment_settings():
    
    varying_params = {}
    
    # --------- Specify parameter values manually -----------------------------
   
    # parameters to keep fixed (but different from standard_value) across combinations
    fixed_params = {
        'infected_individuals_at_start': 100,
        'final_time': 365*3,
        'nucleotide_substitution_rate': 1.105215e-04,
        'long_evo_rate_f': 176.58,
        'R': 1.05,
        'sequencing_rate': 0.1,
        'long_shedders_ratio': 0.01,
        'sequence_long_shedders':True
    }

    n_seeds = 100
    
    return (varying_params,fixed_params,n_seeds)

def main():

    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Run script to test pipeline and OSR values for normal and long shedders")
    parser.add_argument('runner', type=str, help="runner")
    parser.add_argument('experiment_number', type=int, help="experiment number")
    args = parser.parse_args()
    
    # Run the script 
    run_experiment_script(args.runner, 
                          args.experiment_number, 
                          user_set_experiment_settings,
                          experiment_name)
                          
    experiment_numbered_name = f'{experiment_name}_#{args.experiment_number}'
    parameter = 'nucleotide_substitution_rate'
    
    pm.plot_combined_tempest_regressions(experiment_numbered_name, parameter, 
                                      min_seq_number=30, min_sim_lenght=365, individual_type=None,
                                      y_axis_max=0.05)
                                      
    pm.plot_combined_tempest_regressions(experiment_numbered_name, parameter, 
                                      min_seq_number=30, min_sim_lenght=365, individual_type='normal',
                                      y_axis_max=0.05)
                                      
    pm.plot_combined_tempest_regressions(experiment_numbered_name, parameter, 
                                      min_seq_number=0, min_sim_lenght=365, individual_type='long_shedder',
                                      y_axis_max=0.02)

if __name__ == "__main__":
    main()
    
    
