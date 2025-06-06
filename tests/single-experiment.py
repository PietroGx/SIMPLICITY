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
@author: jbescudie
"""

##### <test fixtures>
from tests.common_test_fixtures import fixture_simplicity_runners
from tests.common_test_fixtures import fixture_users_run_seeded_simulation_functions

## fixture  experiment settings (sm.write_settings arguments)
def fixture_experiment_settings():
    parameters      = {'evolutionary rate': [0.001]}
    n_seeds         = 2
    yield (parameters, n_seeds)

##### </test fixtures>

##### <actual test>
def test_run_single_experiment(experiment_name, 
                               experiment_settings, 
                               simplicity_runner, 
                               run_seeded_simulation,
                               plot_trajectory):
    print('')
    print('##########################################')
    print('')

    # import simplicity
    import simplicity.dir_manager           as dm
    import simplicity.settings_manager as sm
    # import simplicity.output_manager   as om
    sr = simplicity_runner

    # set parameters and decide if search or not
    parameters, n_seeds = experiment_settings
    sm.write_settings(parameters, n_seeds)
    # setup experiment files directories
    dm.create_directories(experiment_name)
    # Write experiment settings file
    sm.write_experiment_settings(experiment_name)
    # write simulation parameters files
    sm.read_settings_and_write_simulation_parameters(experiment_name)
    # write seeded simulation parameters files
    sm.write_seeded_simulation_parameters(experiment_name)
    
    # let one of simplicity.runners run each seeded simulation
    sr.run_seeded_simulations(experiment_name, run_seeded_simulation, plot_trajectory)
    
    # archive experiment
    # om.archive_experiment(experiment_name)
    # done
    print('')
    print(f'EXPERIMENT {experiment_name} EXECUTED SUCCESSFULLY.')
    print('')
    print('##########################################')
##### </actual test>


##### run all tests
if __name__ == "__main__":
    
    # all fixtures set combinations
    import itertools
    fixtures_product = itertools.product(
        reversed(tuple(fixture_users_run_seeded_simulation_functions())),
        fixture_experiment_settings(),
        reversed(tuple(fixture_simplicity_runners())),
    )

    # run the test on each fixtures set
    for test_i, fixtures in enumerate(fixtures_product):
        run_seeded_simulation, experiment_settings, simplicity_runner = fixtures
        experiment_name = f'test_{test_i:02d}_single_experiment__{run_seeded_simulation.__name__}__{simplicity_runner.__name__}'        
        test_run_single_experiment(experiment_name, 
                                   experiment_settings, 
                                   simplicity_runner, 
                                   run_seeded_simulation)
