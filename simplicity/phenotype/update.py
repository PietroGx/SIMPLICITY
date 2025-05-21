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

# -*- coding: utf-8 -*-

import simplicity.phenotype.distance as dis
import numpy as np

def immune_waning_fitness_score(population,lineage_genome,consensus):
    # get the hamming distance from the weighted consensus
    distance_from_weighted_consensus = dis.hamming_iw(lineage_genome,consensus)
    # compute fitness score
    infected_fraction = min((population.diagnosed + population.recovered)/population.size,1)
    non_infected_fraction = max((population.size-(population.diagnosed + population.recovered))/population.size,0)
    fitness_score = infected_fraction*distance_from_weighted_consensus + non_infected_fraction/population.active_lineages_n
    
    return fitness_score

def update_fitness_factory(type):
    '''
    Factory of fitness update function. Returns update_fitness, depending on selected
    phenotype model. Update_fitness computes and assigns the  fitness score of every 
    intra host lineage for all individuals to be updated.
    '''
    if type == "linear":
        
        def update_fitness(population,individuals_to_update):
            # individuals - dictionary of individuals in the simulation
            # individuals_to_update - indices of individuals to be updated
            for individual in sorted(individuals_to_update):
                fitness = []
                for lineage_name in population.individuals[individual]['IH_lineages']:
                    lineage_genome = population.get_lineage_genome(lineage_name)
                    fitness.append(dis.hamming(lineage_genome))
                population.individuals[individual]['IH_lineages_fitness_score'] = fitness
                population.individuals[individual]['fitness_score'] = np.average(fitness)
        
        return update_fitness
    
    elif type == "immune_waning":
        
        def update_fitness(population,individuals_to_update,consensus):
            # individuals - dictionary of individuals in the simulation
            # individuals_to_update - indices of individuals to be updated
            # consensus - consensus sequence
            for individual_index in sorted(individuals_to_update):
                fitness = []
                for lineage_name in population.individuals[individual_index]['IH_lineages']:
                    lineage_genome = population.get_lineage_genome(lineage_name)
                    fitness.append(immune_waning_fitness_score(population,lineage_genome,consensus))
                population.individuals[individual_index]['IH_lineages_fitness_score'] = fitness
                population.individuals[individual_index]['fitness_score'] = np.average(fitness)
        
        return update_fitness

