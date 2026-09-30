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

def consensus_distance(lineage_genome, consensus, use_distribution=False):
    # consensus is (sequence, column distributions, delta) from get_consensus
    chat, column_dist, delta = consensus
    if use_distribution:
        return dis.distributional(lineage_genome, column_dist, delta)
    return dis.hamming_iw(lineage_genome, chat)


def immune_waning_fitness_score(population,lineage_genome,consensus,
                                use_distribution=False):
    return fitness_from_distance(
        population,
        consensus_distance(lineage_genome, consensus, use_distribution))


def fitness_from_distance(population, distance_from_weighted_consensus):
    # compute fitness score
    infected_fraction = min((population.diagnosed + population.recovered)/population.size,1)
    non_infected_fraction = max((population.size-(population.diagnosed + population.recovered))/population.size,0)
    
    fitness_score = infected_fraction*distance_from_weighted_consensus + non_infected_fraction/population.active_lineages_n
    
    #epsilon = 1e-6  # fitness floor
    #fitness_score = max(fitness_score, epsilon)
    
    return fitness_score

# def update_relative_fitness(population):
#     infectious_i_list = sorted(population.infectious_i)
   
#     fitness_inf = [population.individuals[i]['fitness_score'] for i in infectious_i_list]
#     fitsum = np.sum(fitness_inf)
    
#     for i in infectious_i_list : 
#         population.individuals[i]['fitness_score'] /= fitsum

def update_fitness_factory(type, consensus_mode='argmax'):
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
                ind = population.individuals[individual]
                scores = {}
                for lineage_name in ind['IH_lineages']:
                    if lineage_name not in scores:
                        scores[lineage_name] = dis.hamming(
                            population.get_lineage_genome(lineage_name))
                ind['IH_lineages_fitness_score'] = [scores[l] for l in ind['IH_lineages']]
                # host score averages over DISTINCT lineages: copies do not count twice
                ind['fitness_score'] = np.average(list(scores.values()))
            # # update relative fitness scores for the population
            # update_relative_fitness(population)
        
        return update_fitness
    
    elif type == "immune_waning":
        use_distribution = consensus_mode == 'distribution'
        # A lineage's distance from the consensus cannot change while that
        # consensus holds: its genome is fixed at creation and a mutation makes
        # a NEW lineage. Keep the distances until the consensus is rebuilt.
        distance_cache = {}
        cached_for = [None]

        def update_fitness(population,individuals_to_update,consensus):
            if cached_for[0] is not consensus:
                distance_cache.clear()
                cached_for[0] = consensus
            # individuals - dictionary of individuals in the simulation
            # individuals_to_update - indices of individuals to be updated
            # consensus - consensus sequence
            for individual_index in sorted(individuals_to_update):
                ind = population.individuals[individual_index]
                scores = {}
                for lineage_name in ind['IH_lineages']:
                    if lineage_name not in scores:
                        d = distance_cache.get(lineage_name)
                        if d is None:
                            d = consensus_distance(
                                population.get_lineage_genome(lineage_name),
                                consensus, use_distribution)
                            distance_cache[lineage_name] = d
                        scores[lineage_name] = fitness_from_distance(population, d)
                ind['IH_lineages_fitness_score'] = [scores[l] for l in ind['IH_lineages']]
                # host score averages over DISTINCT lineages: copies do not count twice
                ind['fitness_score'] = np.average(list(scores.values()))
            # # update relative fitness scores for the population
            # update_relative_fitness(population)
        
        return update_fitness

    raise ValueError(f"Unknown phenotype_model: {type!r}")
