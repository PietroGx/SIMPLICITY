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

def run_seeded_simulation(seeded_simulation_parameters_path: str, 
                          experiment_name: str,
                          plot_trajectory: bool) -> None:
    """runs one seeded simulation. This function isolates how to run one simulation 
    from the looping over all seeded simulation parameters of one experiment.
    
    Support for running Slurm requires this function to be importable from its 
    name as reference (simplicity.runners.unit_run.run_seeded_simulation).
    Also this function's arguments are passed by reference (achieved to 
    serializing/deserializing the reference to the simulation parameters 
    that are resolved via the settings_manager.read_seeded_simulation_parameters function).
    """
    # <settings and output managers>
    import simplicity.settings_manager as sm
    import simplicity.output_manager   as om
    parameters       = sm.read_seeded_simulation_parameters(experiment_name, seeded_simulation_parameters_path)
    output_directory = om.setup_output_directory           (experiment_name, seeded_simulation_parameters_path)
    # create simulation id
    seed = parameters['seed']
    sim_name = sm.generate_filename_from_params(parameters)
    sim_id = f'{experiment_name}: {sim_name}: {seed}'
    
    # </settings and output managers>

    ## <simplicity core>        
    import simplicity.simulation       as sim
    simulation = sim.Simplicity       (parameters, output_directory, sim_id)
    simulation.run()
    ## </simplicity core>        
    if plot_trajectory:    
        simulation.plot()
    