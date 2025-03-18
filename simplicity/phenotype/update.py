# -*- coding: utf-8 -*-

import simplicity.phenotype.distance as dis
import numpy as np

def get_individuals_to_update(subst_coord):
    # get indices of individuals that mutated to update their fitness scores
    mutated_individuals = [i[0] for i in subst_coord]
    mutated_individuals = list(set(mutated_individuals))
    
    return mutated_individuals
 
def update_fitness_factory(type):
    
    if type == "distance from wt":
        
        def update_fitness(population,individuals_to_update):
            # individuals - dictionary of individuals in the simulation
            # individuals_to_update - indices of individuals to be updated
            for individual in individuals_to_update:
                fitness = []
                for lineage_name in population.individuals[individual]['IH_lineages']:
                    lineage_genome = population.get_lineage_genome(lineage_name)
                    fitness.append(dis.hamming(lineage_genome))
                population.individuals[individual]['IH_virus_fitness'] = fitness
                population.individuals[individual]['fitness'] = np.average(fitness)
        
        return update_fitness
    
    elif type == "immune waning":
        
        def update_fitness(population,individuals_to_update,consensus):
            # individuals - dictionary of individuals in the simulation
            # individuals_to_update - indices of individuals to be updated
            # consensus - consensus sequence
            for individual_index in individuals_to_update:
                fitness = []
                for lineage_name in population.individuals[individual_index]['IH_lineages']:
                    lineage_genome = population.get_lineage_genome(lineage_name)
                    fitness.append(dis.hamming_iw(lineage_genome,consensus))
                population.individuals[individual_index]['IH_virus_fitness'] = fitness
                population.individuals[individual_index]['fitness'] = np.average(fitness)
        
        return update_fitness
