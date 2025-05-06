#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 1 11:15:10 2023

@author: pietro
"""
import simplicity.phenotype.update    as pheno
import numpy as np
import copy

# import Evolution.Reference as ref

'''
you can use this code to run the evolution model with variable site substitution rates
you need to pass p_sites_var to the function that selects the substitutions.
You also need to uncomment lines 107,108,109 and pass the p=p_sites_normed 
to the random choice function on line 121

---
avg_spike_ssy = 8.066 * 1e-4

# substitution rates (per site, per year) of SARS-CoV-2 Spike gene
rate_sites_var = [avg_spike_ssy for i in range(len(ref.spike_sequence))]

# substitution rate (per year) of SARS-CoV-2 lineages
E = np.sum(rate_sites_var)

# substitution site relative probability
p_sites_var = [i/E for i in rate_sites_var]
---

'''

def sub_events(rng, e, dt, IH_lineages):
    ''' Choose the number of substitutions that happen in a time step dt,
        drawing random number sampled from Poisson distribution.

    Parameters
    ----------
    rng : numpy method
        Random number generator.
    e : float
        population nucleotide substitution rate.
    dt : float
        delta t - extrande time step.
    IH_lineages : int
        number of active IH_lineages that can mutate in time step dt

    Returns
    -------
    int
        number of substitutions that happen in a time step dt.

    '''
    lambda_pop = e*dt*IH_lineages
    events = rng.poisson(lambda_pop)
    
    return events

def select_positions(population,L_lin,e,dt):
    '''
    Select positions of the pooled genome to mutate.

    Parameters
    ----------
    population : 
        instance of population class.
    L_lin : int
        lenght of IH_lineage genome
    p_sites_var : array
        vector of single site substitution rate
    rng : numpy method
        Random number generator.
    e : float
        population nucleotide substitution rate.
    dt : float
        delta t - extrande time step.

    Raises
    ------
    ValueError
        raise error if the lenght of active_lineages is not same as 
        population.active_lineages_n.

    Returns
    -------
    positions: array
        positions that will undergo substitution
    active_lineages : list
        list of "active" lineages, coord map onto individuals
    sub_number: int
        number of substitution events happening
    '''
    rng = population.rng5
    # list of "active" lineages, coord map onto individuals
    active_lineages = [[i, j] for i in sorted(population.infected_i) for j 
                       in range(0,population.individuals[i]['IH_lineages_number'])]
     
    # raise error if the lenght of active_lineages is not same as active_lineages_n       
    if population.active_lineages_n != len(active_lineages):
        raise ValueError('active_lineages_n must have lenght "active_lineages"')
    
    # vector of pooled genome positions
    positions = np.arange(0,population.active_lineages_n*L_lin,1)
    
    # # vector of pooled genome sites rate
    # p_sites = p_sites_var*population.active_lineages_n
    # p_sites_normed = [i/np.sum(p_sites) for i in p_sites]
    
    # number of substitutions event in time step dt
    sub_number = sub_events(rng, e, dt, population.active_lineages_n)
    
    if sub_number == 0:
        return 'No substitutions'
    else:    
        # return vector of positions that will be mutated p=p_sites_normed
        return rng.choice(positions,sub_number), active_lineages, sub_number

def map_to_lin(positions,L_lin):
    '''
    Maps positions to be mutated onto active_lineages (lineages).

    Parameters
    ----------
    positions: array
        positions that will undergo substitution
    L_lin : int
        lenght of lineage genome

    Returns
    -------
    sub_coord_varmap : tuple
        lineages mutation coordinates:
            (LIST index of active_lineage to mutate,
             position in the genome to mutate).
        the index here does not refer to the index of the population.individuals dic
        

    '''
    sub_coord_linmap = []
    for position in positions:
        c = [int(np.floor(position/L_lin)),position%L_lin]
        sub_coord_linmap.append(c)
        
    return sub_coord_linmap

def map_to_dic(active_lineages, positions, L_lin):
    '''
    Map mutation positions onto dictionary of population.individuals

    Parameters
    ----------
    active_lineages : list
        list of "active" lineages, coord map onto individuals
    positions: array
        positions that will undergo substitution 

    Returns
    -------
    sub_coord_dicmap : list
    
    coordinates of substitution:
        
    index of individual (in population.individuals) that hosts the lineage    
    index of lineage inside the host
    position in the genome of lineage to be mutated
    
    '''
    
    sub_coord_dicmap = []
    sub_coord_linmap = map_to_lin(positions, L_lin)
    
    for i in sub_coord_linmap:
        # i[0] --> index of active_lineages that refers to lineages to mutate 
        lin_coord = active_lineages[i[0]]
        
        # lin_coord[0] --> index of individual (in population.individuals) that hosts the lineage
        # lin_coord[1] --> index of lineage inside the host
        # i[1] --> position in the genome of lineage to be mutated
        sub_coord_dicmap.append([lin_coord[0], 
                           lin_coord[1],
                           i[1]])
    return sub_coord_dicmap
        
def fetch_bases(population, sub_coord_dicmap):
    '''
    Add base to be substituted to subst_coord. First it checks if the position
    is already mutated, if so it takes the mutated base.
    Otherwise fetches the base from the reference genome.

    Parameters
    ----------
    population : instance of population class
    
    sub_coord_dicmap : list
        coordinates of substitution  

    Returns
    -------
    sub_coord_dicmap: list
    
    coordinates of substitution:
        
    index of individual (in population.individuals) that hosts the lineage    
    index of lineage inside the host
    position in the genome of lineage mutated
    nitrogenous base of the lineage at that position

    '''
    # loop over list of coordinates to fetch the corresponding base
    for coord in sub_coord_dicmap:
        # lineage --> list of mutations that a lineage holds (property of individual)        
        lineage_name = population.individuals[coord[0]]['IH_lineages'][coord[1]]
        lineage_genome = population.get_lineage_genome(lineage_name)
        
        # get mutation if the lineage already have one at position
        mut_coord = [mut for mut in lineage_genome if coord[2] in mut]
        # if so, fetch the base from lineage
        if mut_coord:
            coord.append(mut_coord[0][1])
            
        # else, fetch it from reference genome
        else:
            coord.append(population.ref_genome[coord[2]])
        
    return sub_coord_dicmap

def get_subst_numbers(sub_coord_dicmap, sub_number):
    '''
    Returns count of how many of each nitrogenous base will undergo substitution

    Parameters
    ----------
    sub_coord_dicmap : list
        coordinates of substitutions
    sub_number : int
        number of substitution events

    Raises
    ------
    ValueError
        Raise error if the list of coordinates is not the same size as the 
        number of substitution events 

    Returns
    -------
    TYPE
        dic that counts the number of bases to be substituted.
        e.g. 4 A, 3 T, 0 C, 5 G

    '''
    if sub_number != len(sub_coord_dicmap):
        raise ValueError('subst_coord must have sub_number lenght')
    
    A = sum('A' in coord for coord in sub_coord_dicmap)
    T = sum('T' in coord for coord in sub_coord_dicmap)
    C = sum('C' in coord for coord in sub_coord_dicmap)
    G = sum('G' in coord for coord in sub_coord_dicmap)
    # print("bases count", {'A':A,'T':T,'C':C,'G':G})
    return {'A':A,'T':T,'C':C,'G':G} # bases_count

def sub_matrix():
    '''
    Define substitution matrix. If you want to specify a different matrix, change
    the code here.
    '''

    # row A
    pAA = 0 
    pAT = 1/3
    pAC = 1/3
    pAG = 1/3
    
    # row T
    pTA = 1/3
    pTT = 0
    pTC = 1/3
    pTG = 1/3
    
    # row C
    pCA = 1/3
    pCT = 1/3
    pCC = 0
    pCG = 1/3
    
    # row G
    pGA = 1/3
    pGT = 1/3
    pGC = 1/3
    pGG = 0
    
    M = {
        'A': [pAA,pAT,pAC,pAG],
        'T': [pTA,pTT,pTC,pTG],
        'C': [pCA,pCT,pCC,pCG],
        'G': [pGA,pGT,pGC,pGG],
        }
    
    return M
    
def substitution(M,bases_count,rng):
    '''
    Generate substitutions for each base based on substitution matrix. These
    substitutions are not assigned to any genome position yet, the function
    returns a dict with a list of substitutions for each nitrogenous base.

    Parameters
    ----------
    M : substitution matrix
    bases_count : dic
        dictionary with count of bases to substitute

    Returns
    -------
    dic
        dictionary that contain a list of substitutions for each base type.
        e.g A = [T,T,G,C,G,T,T,C,C,G]

    '''
    
    NB = ['A','T','C','G']
    
    A = rng.choice(NB,bases_count['A'],p=M['A'])
    T = rng.choice(NB,bases_count['T'],p=M['T'])
    C = rng.choice(NB,bases_count['C'],p=M['C'])
    G = rng.choice(NB,bases_count['G'],p=M['G'])
    # print("unassigned subst", {'A':A,'T':T,'C':C,'G':G})
    return {'A':A,'T':T,'C':C,'G':G} # unassigned_subst
    
def assign_sub(unassigned_subst,sub_coord_dicmap):
    '''
    Update the substitution coordinates to include the substitutions that will
    be included (assign the substitutions to lineages with the positions in their 
    genome)

    Parameters
    ----------
    unassigned_subst : dic
        dict with a list of substitutions for each nitrogenous base.
    sub_coord_dicmap : coordinates of the substitutions that will be introduced

    Returns
    -------
    sub_coord_dicmap: list
    
    UPDATED coordinates of substitution:
        
    index of individual (in population.individuals) that hosts the lineage    
    index of lineage inside the host
    position in the genome of lineage mutated
    nitrogenous base that will be introduced at that position
    '''
    
    i_A = 0
    i_T = 0
    i_C = 0
    i_G = 0
    
    for coord in sub_coord_dicmap:
        
        if coord[3] == 'A':
            coord[3] = unassigned_subst['A'][i_A]
            i_A +=1
        elif coord[3] == 'T':
            coord[3] = unassigned_subst['T'][i_T]
            i_T +=1
        elif coord[3] == 'C':
            coord[3] = unassigned_subst['C'][i_C]
            i_C +=1

        elif coord[3] == 'G':
            coord[3] = unassigned_subst['G'][i_G]
            i_G +=1
    # print(i_A)
    # print(i_T)
    # print(i_C)
    # print(i_G)    
    return sub_coord_dicmap

def get_lineage_name(phylodots,lineage_name):    
    # function to name the lineages during the simulation
    n_dots = lineage_name.count('.')
    try: 
        phylodots[n_dots] += 1
        lineage_name = lineage_name + '.' + str(phylodots[n_dots])
    except: 
        phylodots.append(0)
        lineage_name = lineage_name + '.0'
    
    return phylodots, lineage_name

def update_lineages(population, sub_coord_dicmap):
    '''
    Update lineages genomes in population with the selected mutations from the 
    coordinates.

    Parameters
    ----------
    population : instance of population class 
    sub_coord_dicmap : substituions coordinates

    '''
    
    # mutate every individual selected and update nodes info
    for coord in sub_coord_dicmap:
        # parent of the phylogenetic tree node
        parent_lineage_name = population.individuals[coord[0]]['IH_lineages'][coord[1]]
        
        # rename the IH lineage after new mutation introduced (in individual)
        population.phylodots, new_lineage_name = get_lineage_name(population.phylodots,parent_lineage_name)
        population.individuals[coord[0]]['IH_lineages'][coord[1]] = new_lineage_name
        # update lineage ih trajectories
        population.individuals[coord[0]]['IH_lineages_trajectory'][new_lineage_name] = {}
        population.individuals[coord[0]]['IH_lineages_trajectory'][new_lineage_name]['ih_birth'] = population.time
        population.individuals[coord[0]]['IH_lineages_trajectory'][new_lineage_name]['ih_death'] = None
        if parent_lineage_name not in population.individuals[coord[0]]['IH_lineages']:
            population.individuals[coord[0]]['IH_lineages_trajectory'][parent_lineage_name]['ih_death'] = population.time
            
    
        # update lineages count inside individual
        population.individuals[coord[0]]['IH_unique_lineages_number'] = len(
            set(population.individuals[coord[0]]['IH_lineages']))
        
        parent_lineage_genome = population.get_lineage_genome(parent_lineage_name)
        new_lineage_genome = copy.deepcopy(parent_lineage_genome)
        new_lineage_genome.append([coord[2],coord[3]])
        new_lineage_genome.sort(key=lambda x: x[0])  # sort by position

        population.phylogenetic_data.append(
                     {'Time_emergence'  : population.time,
                      'Lineage_name'    : new_lineage_name,
                      'Lineage_parent'  : parent_lineage_name,
                      'Genome'          : new_lineage_genome,
                      'Host_type'       : population.individuals[coord[0]]['type']
                    })
        # print(f'Phylo data: {population.phylogenetic_data}')
    mutated_ids = {coord[0] for coord in sub_coord_dicmap}

    for idx in mutated_ids:
        indiv = population.individuals[idx]
        combined = sorted(zip(indiv['IH_lineages'], indiv['IH_lineages_fitness_score']))
        indiv['IH_lineages'], indiv['IH_lineages_fitness_score'] = map(list, zip(*combined))
        
def get_individuals_to_update(subst_coord):
    # get indices of individuals that mutated to update their fitness scores
    mutated_individuals = [i[0] for i in subst_coord]
    mutated_individuals = sorted(set(mutated_individuals))
    
    return mutated_individuals

def mutate(population, e, dt, phenotype_model, *args):
    '''
    Mutation model, mutates the viruses in the population.

    Parameters
    ----------
    e : float
        population nucleotide substitution rate.
    dt : float
        delta t - extrande time step (in years)

    '''
    L_lin = population.L_lin
    # select number of substitutions and positions in pooled genome
    select_pos = select_positions(population, L_lin, e, dt) 
    
    # if mutations are happening
    if select_pos != 'No substitutions':
        
        # positions in pooled genome
        positions = select_pos[0] 
        # vector of active lineages
        active_lineages = select_pos[1]
        # number of substitutions happening
        subst_number = select_pos[2]
        
        # map substitutions to index of lineages
        subst_coord = map_to_dic(active_lineages, positions, L_lin)
        
        # fetch the bases that will undergo substituion
        subst_coord = fetch_bases(population, subst_coord)
        # count number of each nitrogenous base to be mutated (for bulk update)
        bases_count = get_subst_numbers(subst_coord, subst_number)
        # use substitution matrix to select mutations that will happen
        unassigned_subst = substitution(sub_matrix(), bases_count,population.rng5)
        # assign the mutations to their relative genome
        subst_coord = assign_sub(unassigned_subst, subst_coord)
        # update lineages in simulation with the corresponding substitution
        update_lineages(population, subst_coord)
        
        # update fitness score of individuals and lineages
        
        # print('--X----X----X----X----X----X--')
        # print('Substitution coordinates: ', subst_coord)
        individuals_to_update = get_individuals_to_update(subst_coord)
        # print('Individuals to be updated: ', individuals_to_update)
        update_fitness = pheno.update_fitness_factory(phenotype_model)
        if args:
            consensus = args[0]
            update_fitness(population, individuals_to_update, consensus)
        else:
            update_fitness(population, individuals_to_update)

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        