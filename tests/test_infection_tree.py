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


## fixture  experiment settings (sm.write_settings arguments)
def fixture_experiment_settings():
    parameters      = {'population_size': [20],"infected_individuals_at_start": [3]}
    n_seeds         = 1
    return (parameters, n_seeds)

def test_run_single_seeded_simulation(experiment_name, 
                                      experiment_settings):
    
    # import simplicity
    import simplicity.dir_manager           as dm
    import simplicity.settings_manager as sm
    import simplicity.runners.unit_run as run
    # setup experiment files directories
    dm.create_directories(experiment_name)
    
    # set parameters 
    parameters, n_seeds = experiment_settings
    sm.write_settings(parameters, n_seeds)
    # Write experiment settings file
    sm.write_experiment_settings(experiment_name)
    # write simulation parameters files
    sm.read_settings_and_write_simulation_parameters(experiment_name)
    # write seeded simulation parameters files
    sm.write_seeded_simulation_parameters(experiment_name)
    
    # get seeded simulation parameters paths and take only the first entry
    seeded_simulation_parameters_path = sm.get_seeded_simulation_parameters_paths(experiment_name)[0]
    
    # run single seeded simulation 
    run.run_seeded_simulation(seeded_simulation_parameters_path, experiment_name)
    

    print('')
    print(f'EXPERIMENT {experiment_name} EXECUTED SUCCESSFULLY.')
    print('')
    print('##########################################')


if __name__ == "__main__":
    
    test_run_single_seeded_simulation('Test infection tree #14', 
                                      fixture_experiment_settings())