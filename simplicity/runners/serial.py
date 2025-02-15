#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: jbescudie

"""
import simplicity.dir_manager as dm


def run_seeded_simulations(experiment_name, run_seeded_simulation, plot_trajectory):
    seeded_simulation_parameters_paths = dm.get_seeded_simulation_parameters_paths(experiment_name)
    # run a Simplicity simulation for each seeded parameters
    for seeded_simulation_parameters_path in seeded_simulation_parameters_paths:
        run_seeded_simulation(seeded_simulation_parameters_path, experiment_name, plot_trajectory)
