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
import argparse
import numpy as np
experiment_name =  'OSR_fit'

def fixture_experiment_settings():
    # number of values for NSR
    nucleotide_substitution_rate_num_values = 15
    
    # Generate values spaced logarithmically between 10^-5 and 10^-2
    values = np.logspace(np.log10(1e-6), np.log10(0.0003), 
                         num=nucleotide_substitution_rate_num_values)
    nucleotide_substitution_rate_values = values.tolist()
    # parameters value to get combinations from
    varying_params = {
        'nucleotide_substitution_rate': nucleotide_substitution_rate_values
    }
    # parameters to keep fixed (but different from standard_value) across combinations
    fixed_params = {
        'infected_individuals_at_start': 10,
        'final_time': 365*3,
        'R': 1.1
    }
    
    n_seeds = 100
    
    return (varying_params,fixed_params,n_seeds)

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Run script to generate IH lineages data")
    parser.add_argument('runner', type=str, help="runner")
    parser.add_argument('experiment_number', type=int, help="experiment number")
    args = parser.parse_args()
    # Run the script 
    run_experiment_script(args.runner, 
                          args.experiment_number, 
                          fixture_experiment_settings,
                          experiment_name)

if __name__ == "__main__":
    main()
    
    
