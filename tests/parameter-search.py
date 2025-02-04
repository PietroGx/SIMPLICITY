#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: jbescudie
"""

##### <test fixtures>
from simplicity.tests.common_test_fixtures import fixture_simplicity_runners
from simplicity.tests.common_test_fixtures import fixture_users_run_seeded_simulation_functions

## fixture  experiment settings (sm.write_settings arguments)
def fixture_experiment_settings():
    parameters      = {'evolutionary rate': [0.001]}
    n_seeds         = 10
    # parameter search
    target_value    = 0.002  
    yield (parameters, n_seeds, target_value)

def fixture_parameter_search_settings():
    THRESHOLD       = 0.0001
    yield (THRESHOLD,)

##### </test fixtures>

##### <actual test>
def test_run_parameter_search(experiment_name, 
                              experiment_settings, 
                              paremeter_search_settings, 
                              simplicity_runner, 
                              run_seeded_simulation):
    # import simplicity
    import simplicity.config                     as config
    import simplicity.settings_manager           as sm
    # import simplicity.output_manager             as om
    import simplicity.tuning.parameter_search    as ps
    #import simplicity.runners.slurm   as sr      # is passed as a fixture instead
    sr = simplicity_runner

    # set parameters 
    parameters, n_seeds, target_value = experiment_settings
    sm.write_settings(parameters, n_seeds, target_value)
    
    # setup experiment files directories
    config.create_directories(experiment_name)
    
    # Write experiment settings file
    sm.write_experiment_settings(experiment_name)
    # write simulation parameters files
    sm.read_settings_and_write_simulation_parameters(experiment_name)
    
    # write parameter search settings
    threshold, = paremeter_search_settings
    ps.write_paremeter_search_settings(experiment_name, threshold)
    # parameter search fit
    ps.fit(experiment_name, sr, run_seeded_simulation)
    
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
    # must be set before first import :/
    import os
    SIMPLICITY_DATA_DIR = os.environ["SIMPLICITY_DATA_DIR"] = os.path.join(os.getcwd(), "Data-tests")

    # all fixtures set combinations
    import itertools
    fixtures_product = itertools.product(
        reversed(tuple(fixture_users_run_seeded_simulation_functions())),
        fixture_experiment_settings(),
        reversed(tuple(fixture_simplicity_runners())),
        fixture_parameter_search_settings(),
    )

    # run the test on each fixtures set
    for test_i, fixtures in enumerate(fixtures_product):
        run_seeded_simulation, experiment_settings, simplicity_runner, parameter_search_settings = fixtures
        experiment_name = f'test_{test_i:02d}_parameter_search__{run_seeded_simulation.__name__}__{simplicity_runner.__name__}'
        test_run_parameter_search(SIMPLICITY_DATA_DIR, experiment_name, experiment_settings, parameter_search_settings, simplicity_runner, run_seeded_simulation)
