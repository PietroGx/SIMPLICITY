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
Created on Fri Mar 28 12:51:33 2025

@author: pietro
"""
import numpy as np

def SIDR_propensities(population, beta_normal, beta_long, k_d, k_v, seq_rate):
    propensities_params = [beta_normal, beta_long, k_d, k_v]
    # system propensities
    # reaction_id, reaction_rate, action
    propensities = [
        ('infection_normal', beta_normal * population.infectious_normal, lambda: infection(population, from_long_shedder=False)),
        ('infection_long',   beta_long * population.infectious_long  ,   lambda: infection(population, from_long_shedder=True)),
        ('diagnosis',        k_d * population.detectables,               lambda: diagnosis(population, seq_rate)),
        ('add_ih_lineage',   k_v * population.infected,                  lambda: add_lineage(population))
    ]
    # print('a1:',beta * population.infectious)
    # print('a2:',k_d * population.detectables)
    # print('a3:',k_v * population.infected)
    return propensities, propensities_params

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
        # patient_id, time, genome, subst number, patient type, infection duration, ih lineage number
        i = 0
        for lineage_name in population.individuals[diagnosed_individual_i]['IH_lineages']:
            genome = population.get_lineage_genome(lineage_name)
            population.sequencing_data.append({
                'individual_index': diagnosed_individual_i,
                'sequencing_time' : population.time,
                'lineage_name'    : lineage_name,
                'sequence'        : genome,
                'sequence_lenght' : len(genome),
                'individual_type' : population.individuals[diagnosed_individual_i]['type'],
                'infection_duration' : population.time - population.individuals[diagnosed_individual_i]['t_infection'],
                'intra-host_lineage_index' :i
            })
            i += 1

    # update the active lineages number
    population.active_lineages_n -= population.individuals[diagnosed_individual_i]['IH_lineages_number']
    
    population.diagnosed += 1
    # set patient as diagnosed
    population.individuals[diagnosed_individual_i]['t_not_infectious'] = population.time
    population.individuals[diagnosed_individual_i]['state'] = 'diagnosed'

   
    if population.individuals[diagnosed_individual_i]['t_not_infected'] is None:
        population.individuals[diagnosed_individual_i]['t_not_infected'] = population.time
    else:
        raise ValueError('Individual already recovered!!')
    
    # remove individual index from detectable_i, infectious_i, infected_i and add to diagnosed_i
    
    if diagnosed_individual_i in population.infectious_i:
        population.infectious_i.remove(diagnosed_individual_i)
        population.infectious -= 1
    
    population.detectable_i.remove(diagnosed_individual_i)
    population.detectables -= 1
    
    # update diagnosed and infected
    population.infected_i.discard(diagnosed_individual_i)
    population.diagnosed_i.add(diagnosed_individual_i)

    # add a susceptible back in from reservoir
    population.susceptibles += 1
    new_susceptible_index = population.reservoir_i.pop()
    population.susceptibles_i.add(new_susceptible_index)

def infect_long_shedder(population, new_infected_index):
    
    min_time = 0
    min_infected_n = 0
    max_long_shedders = population.size * population.long_shedders_ratio
    
    if (population.time > min_time and 
        population.infected > min_infected_n and
        population.long_shedders < max_long_shedders and 
        population.rng4.uniform() < 0.01
                      ):
        population.long_shedder_i.add(new_infected_index)
        population.long_shedders += 1
        
        individual_type = 'long_shedder'
        
    
    else:
        individual_type = 'normal'
        
    return individual_type

def infection(population, from_long_shedder=False):
    '''
    Select a random susceptible individual to be infected and tags it as 
    such. Update compartments and infected_i
    '''
    # Convert sets to lists for random sampling
    if from_long_shedder:
        # Infectious AND long-shedder
        infectious_i_list = sorted(population.infectious_i & population.long_shedder_i)
        parent = population.rng4.choice(infectious_i_list)
    else: 
        # Infectious but NOT long-shedder
        infectious_i_list = sorted(population.infectious_i - population.long_shedder_i)
        infection_fitness = [population.individuals[i]['fitness_score'] for i in infectious_i_list]
        fitness_sum = sum(infection_fitness)
        
        if fitness_sum > 0:
            weights = [f / fitness_sum for f in infection_fitness]
        else:
            # fallback to uniform probabilities if all fitness scores are zero
            weights = [1 / len(infectious_i_list)] * len(infectious_i_list)
        
        parent = population.rng4.choice(infectious_i_list, p=weights)
        
    susceptibles_list = sorted(population.susceptibles_i)
    
    # select random patient to be infected
    new_infected_index = population.rng4.choice(susceptibles_list)
    population.exclude_i = {new_infected_index}
    new_inf = population.individuals[new_infected_index]

    # Move individual from susceptibles to infected
    population.susceptibles_i.discard(new_infected_index)
    population.infected_i.add(new_infected_index)

    # time of infection
    new_inf['t_infection'] = population.time
    # state
    new_inf['state'] = 'infected'
    
    if population.update_ih_mode == 'jump':
        # Sample first jump time from exponential distribution
        individual_type = new_inf['type']
        rate = - population.host_model[individual_type].A[0][0]
        new_inf['t_next_state'] = (
        population.time + population.rng3.exponential(scale=1 / rate)
        )
        
    # parent
    new_inf['parent'] = parent

    # Select random lineage from parent
    index = population.rng4.integers(0, population.individuals[parent]['IH_lineages_number'])
    transmitted_lineage = population.individuals[parent]['IH_lineages'][index]
    transmitted_fitness = population.individuals[parent]['IH_lineages_fitness_score'][index]
    # count lineage infection
    population._phylo_name_map[transmitted_lineage]['Total_infections'] = population._phylo_name_map[transmitted_lineage].get('Total_infections', 0) + 1
    
    # assign transmitted lineage
    new_inf['inherited_lineage']  = transmitted_lineage
    new_inf['IH_lineages_number'] += 1
    # update lineage trajectory
    new_inf['IH_lineages_trajectory'][transmitted_lineage] = {'ih_birth':None,'ih_death':None}
    # update mutation weight event timer
    new_inf['time_last_weight_event'] = population.time
    
    # update the active lineages number
    population.active_lineages_n += new_inf['IH_lineages_number']
    
    # Append transmitted lineage + fitness
    new_lineages = new_inf['IH_lineages']
    new_fitness = new_inf['IH_lineages_fitness_score']
    new_lineages.append(transmitted_lineage)
    new_fitness.append(transmitted_fitness)

    # Sort lineage and fitness together
    combined = sorted(zip(new_lineages, new_fitness))
    lineages_sorted, fitness_sorted = zip(*combined)
    new_inf['IH_lineages'] = list(lineages_sorted)
    new_inf['IH_lineages_fitness_score'] = list(fitness_sorted)
    
    # Update individual fitness (average of lineage's fitness)
    new_inf['fitness_score'] = round(np.average(fitness_sorted), 4)

    # store infection info for R effective
    population.individuals[parent]['new_infections'].append({
        'time_infection': population.time,
        'transmitted_lineage': transmitted_lineage,
        'individual_infected': new_infected_index
    })
    
    # set the individual to be a long shedder if conditions apply
    new_inf['type'] = infect_long_shedder(population, new_infected_index)
    if new_inf['type']  == 'long_shedder':
        new_inf['IH_lineages_max'] = population.rng3.integers(5,16)
    
    # update susceptibles and infected 
    population.susceptibles -= 1
    population.infected += 1

def add_lineage(population):
    """Adds a new lineage to a randomly selected infected individual, 
    duplicating an existing lineage or replacing one if at max capacity."""
    # Convert set to list for sampling
    infected_list = sorted(population.infected_i)
    individual_index = population.rng4.choice(infected_list)

    individual = population.individuals[individual_index]
    IH_lineage_n = individual['IH_lineages_number']
    IH_lineages_max = individual['IH_lineages_max']

    # Add a duplicate lineage if under max capacity 
    if IH_lineage_n < IH_lineages_max:
        idx = population.rng4.integers(0, IH_lineage_n)
        duplicated_lineage = individual['IH_lineages'][idx]
        individual['IH_lineages'].append(duplicated_lineage)
        individual['IH_lineages_fitness_score'].append(individual['IH_lineages_fitness_score'][idx])
        individual['IH_lineages_number'] += 1
        population.active_lineages_n += 1
        
    # Replace one lineage (delete + duplicate) if at max
    elif IH_lineage_n == IH_lineages_max and IH_lineages_max > 1:
        # Randomly delete one
        delete_idx = population.rng4.integers(0, IH_lineage_n)
        deleted_lineage = individual['IH_lineages'].pop(delete_idx)
        individual['IH_lineages_fitness_score'].pop(delete_idx)
        # update lineage trajectory
        if deleted_lineage not in individual['IH_lineages']:
            individual['IH_lineages_trajectory'][deleted_lineage]['ih_death'] = population.time
        # Duplicate another
        idx = population.rng4.integers(0, IH_lineage_n - 1)  # now virus_n - 1 after deletion
        individual['IH_lineages'].append(individual['IH_lineages'][idx])
        individual['IH_lineages_fitness_score'].append(individual['IH_lineages_fitness_score'][idx])
        

    # Keep IH_lineages and IH_lineages_fitness_score sorted
    combined = sorted(zip(individual['IH_lineages'], individual['IH_lineages_fitness_score']))
    individual['IH_lineages'], individual['IH_lineages_fitness_score'] = map(list, zip(*combined))

    # update individual fitness
    individual['fitness_score'] = round(np.average(individual['IH_lineages_fitness_score']), 4)
