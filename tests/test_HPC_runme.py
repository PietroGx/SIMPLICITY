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
Created on Sat Jan 25 11:07:15 2025

@author: pietro
"""
from simplicity.runme import run_experiment
from tests.test_local_runme import test_experiment_output
import simplicity.runners.slurm
import argparse

## fixture  experiment settings (sm.write_settings arguments)
def fixture_experiment_settings():
    parameters      = {}
    n_seeds         = 10
    return (parameters, n_seeds)

##### <actual test>
def test_run_experiment_HPC(test_number:int):

    print('')
    print('##########################################')
    print('testing HPC runner')
    print('##########################################')
    experiment_name = f'test_HPC_experiment_#{test_number}'
    try:
        run_experiment(experiment_name, 
                       fixture_experiment_settings,             
                       simplicity_runner  = simplicity.runners.slurm,
                       plot_trajectory = True,
                       archive_experiment = False)
    except:
        raise RuntimeError('The code did not pass the running test')
    print('')
    print(f'TEST HPC RUNNER #{test_number} -- SUCCESS.')
    print('##########################################')
    return experiment_name
    
def test_experiment_HPC(test_number:int):
    test_experiment_output(test_run_experiment_HPC(test_number))
    
##### </actual test>

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Run test experiment on HPC")
    parser.add_argument('param', type=int, help="Test number")
    args = parser.parse_args()
    # Run the test with the provided parameter
    test_experiment_HPC(args.param)

if __name__ == "__main__":
    main()