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

