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

## fixture simplicity.runners implementation
def fixture_simplicity_runners():
    # serial
    import simplicity.runners.serial          as sr_serial
    yield sr_serial

    # set environ variables (not used by sr_serial)
    import os
    SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS = os.environ["SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS"] = str(5)
    
    # multiprocessing
    import simplicity.runners.multiprocessing as sr_mp
    yield sr_mp
    
    # # slurm
    # import simplicity.runners.slurm           as sr_slurm
    # yield sr_slurm

## fixture user's run_seeded_simulation function (see simplicity.runners)
def fixture_users_run_seeded_simulation_functions():
    # ! note that we respect the import pattern required by simplicity.runners.slurm
    
    # included in simplicity
    from simplicity.runners.unit_run import run_seeded_simulation
    yield run_seeded_simulation
    
    # user provided
    from  tests.common_test_fixtures import user_run_seeded_simulation
    yield user_run_seeded_simulation

## user's run_seeded_simulation function needs to be accessible by getattr(module, fn_name) for simplicity.runners.slurm
def user_run_seeded_simulation(seeded_simulation_parameters_path: str, 
                               experiment_name: str,
                               plot_trajectory: bool) -> None:
    print("Run seeded simulation routine:", seeded_simulation_parameters_path)
    
    # import simplicity

    # <settings and output managers>
    import simplicity.settings_manager as sm
    import simplicity.output_manager   as om
    parameters       = sm.read_seeded_simulation_parameters(experiment_name, seeded_simulation_parameters_path)
    output_directory = om.setup_output_directory           (experiment_name, seeded_simulation_parameters_path)
    # </settings and output managers>

    ## <simplicity core>        
    import simplicity.simulation       as sim
    import simplicity.population       as pop
    simulation = sim.Simplicity       (parameters, output_directory)
    population = pop.create_population(parameters)
    simulation.run(population)
    if plot_trajectory:    
        simulation.plot()
    
    ## </simplicity core>        

    print("Completed seeded simulation routine:", seeded_simulation_parameters_path)

def fixture_tree_builder_type():
    yield 'infection'
    yield 'phylogenetic'
    
def fixture_infection_tree_subtype():
    yield 'binary'
    yield 'compact'
    yield 'fitness'
    
def fixture_phylogenetic_tree_subtype():
    yield 'binary'
    yield 'non-binary'
    
    
    
    
    
    
    
    
    
    
    
    