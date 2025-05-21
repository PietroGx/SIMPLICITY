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
@author: jbescudie
"""
import simplicity.extrande               as e
import simplicity.population             as pop
from   simplicity.random_gen             import randomgen
import simplicity.output_manager         as om
import simplicity.plots_manager          as pm

class Simplicity:
    def __init__(self, parameters, output_directory, sim_id):
        self.parameters       = parameters
        self.output_directory = output_directory
        self.population       = pop.create_population(parameters)
        self.sim_id  = sim_id
       
    def run(self):
        assert not hasattr(self, "simulation_output"), "cannot run simulation twice"
        
        # create random number generators
        seeds_generator=randomgen(self.parameters["seed"])
        rng1 = randomgen(seeds_generator.integers(0,10000)) # for drawing time to next reaction
        rng2 = randomgen(seeds_generator.integers(0,10000)) # for rejection sampling of reactions
        
        # use factory to assign simulation algorithm
        self.extrande = e.extrande_factory(self.parameters["phenotype_model"], 
                                           self.parameters, self.sim_id,
                                           rng1, rng2)
        
        # run simulation
        self.simulation_output = self.extrande(self.population)
        
        # simulation data output
        # print('Saving simulation trajectory data...')
        om.save_simulation_trajectory(self.simulation_output, 
                                      self.output_directory)
        # print('DONE.')
        # print('')
        # print('Saving lineage frequency data...')
        om.save_lineage_frequency(self.simulation_output, 
                                      self.output_directory)
        # print('DONE.')
        # print('')
        # print('Saving sequencing dataset...')
        om.save_sequencing_dataset(self.simulation_output, 
                                      self.output_directory)
        # print('DONE.')
        # print('')
        # print('Saving lineage individuals data...')
        om.save_individuals_data(self.simulation_output, 
                                      self.output_directory)
        # print('DONE.')
        # print('')
        # print('Saving phylogenetic data...')
        om.save_phylogenetic_data(self.simulation_output, 
                                      self.output_directory)
        # print('DONE.')
        # print('')
        # print('Saving fitness trajectory data...')
        om.save_fitness_trajectory(self.simulation_output, 
                                      self.output_directory)
        # print('DONE.')
        # print('')
        # print('Saving final time datapoint...')
        om.save_final_time(self.simulation_output, 
                                      self.output_directory)
        # print('DONE.')
        # print('')
        # om.save_DEBUG_update_ih(self.simulation_output, 
        #                               self.output_directory)
        
    def plot(self):
        assert hasattr(self, "simulation_output"), "simulation has not run yet. Hint: see Simplicity.run"
        # plot simulation trajectory
        pm.plot_simulation(self.output_directory,threshold=0.1)
    
    def save_consensus(self):
        assert hasattr(self, "simulation_output"), "simulation has not run yet. Hint: see Simplicity.run"
        # save consensus sequences calculated during simulation
        # function in om need to be written if I implement this functionality
        # takes population.consensus_sequences_t and writes it to csv
        # also write function to read it into df!
        om.save_consensus_sequences_t(self.simulation_output,
                                      self.output_directory)
        
        
        