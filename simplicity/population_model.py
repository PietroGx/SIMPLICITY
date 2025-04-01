#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 28 12:51:33 2025

@author: pietro
"""
import numpy as np

def SIDR_propensities(population, beta, k_d, k_v, seq_rate):
    propensities_params = [beta, k_d, k_v]
    return [
        (beta * population.infectious_normal, lambda: infection(population)),
        (k_d * population.detectables, lambda: diagnosis(population, seq_rate)),
        (k_v * population.infected, lambda: add_variant(population))
    ], propensities_params

def diagnosis(population, seq_rate=0):
    '''
    Select an infected (and detectable) individual at random and tags it
    as "diagnosed". Update the infected and diagnosed compartments. 
    Update detectable_i, infectious_i, infected_i and diagnosed_i.
    '''
    population.infected -= 1
    
    # select random patient to be diagnosed
    diagnosed_individual_i = population.rng4.choice(sorted(population.detectable_i))
    
    if population.rng6.uniform(0, 1) < seq_rate:
        # store sequencing data
        # patient_id, time, genome, subst number, patient type, infection duration
        i = 0
        for lineage_name in population.individuals[diagnosed_individual_i]['IH_lineages']:
            genome = population.get_lineage_genome(lineage_name)
            population.sequencing_data.append([
                diagnosed_individual_i,
                population.time,
                genome,
                len(genome),
                population.individuals[diagnosed_individual_i]['type'],
                population.time - population.individuals[diagnosed_individual_i]['t_infection'],
                i
            ])
            i += 1

    # update the active variants number
    population.active_variants_n -= population.individuals[diagnosed_individual_i]['IH_virus_number']
    
    population.diagnosed += 1
    # set patient as diagnosed
    population.individuals[diagnosed_individual_i]['t_not_infectious'] = population.time
    population.individuals[diagnosed_individual_i]['state'] = 'diagnosed'
    
    # remove individual index from detectable_i, infectious_i, infected_i and add to diagnosed_i
    if diagnosed_individual_i in population.infectious_normal_i:
        population.infectious_normal_i.remove(diagnosed_individual_i)
        population.infectious_normal -= 1
    
    population.detectable_i.remove(diagnosed_individual_i)
    population.detectables -= 1
    # update diagnosed and infected
    population.infected_i.discard(diagnosed_individual_i)
    population.diagnosed_i.add(diagnosed_individual_i)

    # add a susceptible back in from reservoir
    population.susceptibles += 1
    new_susceptible_index = population.reservoir_i.pop()
    population.susceptibles_i.add(new_susceptible_index)


def infection(population):
    '''
    Select a random susceptible individual to be infected and tags it as 
    such. Update compartments and infected_i
    '''
    # Convert sets to lists for random sampling
    infectious_i_list = sorted(population.infectious_normal_i)
    susceptibles_list = sorted(population.susceptibles_i)

    # select random patient to be the transmitter
    fitness_inf = [population.individuals[i]['fitness'] for i in infectious_i_list]
    fitsum = np.sum(fitness_inf)

    if fitsum > 0:
        normed_fitness = [f / fitsum for f in fitness_inf]
        parent = population.rng4.choice(infectious_i_list, p=normed_fitness)
    else:
        parent = population.rng4.choice(infectious_i_list)

    # select random patient to be infected
    new_infected_index = population.rng4.choice(susceptibles_list)
    population.exclude_i = {new_infected_index}

    # update the active variants number
    population.active_variants_n += population.individuals[new_infected_index]['IH_virus_number']

    # Move individual from susceptibles to infected
    population.susceptibles_i.discard(new_infected_index)
    population.infected_i.add(new_infected_index)

    # time of infection
    population.individuals[new_infected_index]['t_infection'] = population.time
    # state
    population.individuals[new_infected_index]['state'] = 'infected'
    # Sample first jump time from exponential distribution
    rate = - population.individuals[new_infected_index]['model'].A[0][0]
    population.individuals[new_infected_index]['t_next_state'] = (
    population.time + population.rng3.exponential(scale=1 / rate)
    )
    # parent
    population.individuals[new_infected_index]['parent'] = parent

    # Select random lineage from parent
    index = population.rng4.integers(0, population.individuals[parent]['IH_virus_number'])
    transmitted_lineage = population.individuals[parent]['IH_lineages'][index]
    transmitted_fitness = population.individuals[parent]['IH_virus_fitness'][index]

    # Append transmitted lineage + fitness
    new_lineages = population.individuals[new_infected_index]['IH_lineages']
    new_fitness = population.individuals[new_infected_index]['IH_virus_fitness']
    new_lineages.append(transmitted_lineage)
    new_fitness.append(transmitted_fitness)

    # Sort lineage and fitness together
    combined = sorted(zip(new_lineages, new_fitness))
    lineages_sorted, fitness_sorted = zip(*combined)
    population.individuals[new_infected_index]['IH_lineages'] = list(lineages_sorted)
    population.individuals[new_infected_index]['IH_virus_fitness'] = list(fitness_sorted)

    # Update individual fitness (average of variants' fitness)
    population.individuals[new_infected_index]['fitness'] = round(np.average(fitness_sorted), 4)

    # store infection info for R effective
    population.individuals[parent]['new_infections'].append({
        'time_infection': population.time,
        'transmitted_lineage': transmitted_lineage
    })

    # update susceptibles and infected 
    population.susceptibles -= 1
    population.infected += 1

def add_variant(population):
    """Adds a new variant to a randomly selected infected individual, 
    duplicating an existing lineage or replacing one if at max capacity."""
    # Convert set to list for sampling
    infected_list = sorted(population.infected_i)
    individual_index = population.rng4.choice(infected_list)

    individual = population.individuals[individual_index]
    virus_n = individual['IH_virus_number']
    virus_max = individual['IH_virus_max']

    # Add a duplicate lineage if under max capacity 
    if virus_n < virus_max:
        idx = population.rng4.integers(0, virus_n)
        individual['IH_lineages'].append(individual['IH_lineages'][idx])
        individual['IH_virus_fitness'].append(individual['IH_virus_fitness'][idx])
        individual['IH_virus_number'] += 1
        population.active_variants_n += 1
    # Replace one variant (delete + duplicate) if at max
    elif virus_n == virus_max and virus_max > 1:
        # Randomly delete one
        delete_idx = population.rng4.integers(0, virus_n)
        individual['IH_lineages'].pop(delete_idx)
        individual['IH_virus_fitness'].pop(delete_idx)
        # Duplicate another
        idx = population.rng4.integers(0, virus_n - 1)  # now virus_n - 1 after deletion
        individual['IH_lineages'].append(individual['IH_lineages'][idx])
        individual['IH_virus_fitness'].append(individual['IH_virus_fitness'][idx])

    # Keep IH_lineages and IH_virus_fitness sorted
    combined = sorted(zip(individual['IH_lineages'], individual['IH_virus_fitness']))
    individual['IH_lineages'], individual['IH_virus_fitness'] = map(list, zip(*combined))

    # update individual fitness
    individual['fitness'] = round(np.average(individual['IH_virus_fitness']), 4)
