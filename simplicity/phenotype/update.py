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
        
        def update_fitness(individuals,individuals_to_update, consensus):
            # individuals - dictionary of individuals in the simulation
            # individuals_to_update - indices of individuals to be updated
            # consensus - consensus sequence
            for individual in individuals_to_update:
                
                fitness = []
                for lineage in individuals[individual]['viral_genomes']:
                    fitness.append(dis.hamming(lineage))
                individuals[individual]['IH_virus_fitness'] = fitness
                individuals[individual]['fitness'] = np.average(fitness)
            return individuals
        
        return update_fitness
    
    elif type == "immune waning":
        
        def update_fitness(individuals,individuals_to_update,consensus):
            for individual in individuals_to_update:
                
                fitness = []
                for lineage in individuals[individual]['viral_genomes']:
                    fitness.append(dis.hamming_iw(lineage,consensus))
                individuals[individual]['IH_virus_fitness'] = fitness
                individuals[individual]['fitness'] = np.average(fitness)
            return individuals
        
        return update_fitness
