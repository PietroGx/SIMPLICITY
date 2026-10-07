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

def run_seeded_simulation(experiment_name: str, group: str, index: int) -> None:
    """Runs one repeat. This isolates how to run one repeat from the looping
    over all of an experiment's repeats.

    Support for running Slurm requires this function to be importable from its
    name as reference (simplicity.runners.unit_run.run_seeded_simulation).
    Its arguments are passed by value through the environment, which is why a
    repeat is now (group, index) rather than a path: the task resolves both
    through simplicity.jobs instead of reconstructing a filesystem location.
    """
    # <settings and output managers>
    import simplicity.jobs            as jobs
    import simplicity.settings_manager as sm
    import simplicity.output_manager   as om

    repeat = jobs.get_repeat(experiment_name, group, index)
    # the repeat's only contribution to its own parameters is the seed
    parameters = sm.read_simulation_parameters(experiment_name, repeat['stem'])
    parameters['seed'] = repeat['seed']

    output_directory = om.setup_output_directory(experiment_name, group, repeat)
    # create simulation id. The stem, not generate_filename_from_params over the
    # complete parameter set -- that printed ~250 characters of every parameter
    # on every line, and the stem already names the simulation and its id.
    sim_id = f'{experiment_name}: {group}/{repeat["stem"]}: {repeat["seed"]}'

    # </settings and output managers>

    ## <simplicity core>
    import simplicity.simulation       as sim
    # Its own file, not the state record: a running simulation rewrites this
    # every 30s, and a read-modify-write of the state could land after a
    # terminal COMPLETED and put STARTED back.
    progress_file_path = jobs.progress_path(experiment_name, group, index)
    simulation = sim.Simplicity       (parameters, output_directory, sim_id,
                                       progress_file_path=progress_file_path)
    simulation.run()
    ## </simplicity core>
