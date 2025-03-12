#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 1 11:15:10 2023

@author: pietro
"""
import simplicity.phenotype.distance as dis
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

# substitution rate (per year) of SARS-CoV-2 variants
E = np.sum(rate_sites_var)

# substitution site relative probability
p_sites_var = [i/E for i in rate_sites_var]
---

'''

def sub_events(rng, e, dt, variants):
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
    variants : int
        number of active variants that can mutate in time step dt

    Returns
    -------
    int
        number of substitutions that happen in a time step dt.

    '''
    lambda_pop = e*dt*variants
    events = rng.poisson(lambda_pop)
    
    return events

def select_positions(population,L_var,rng,e,dt):
    '''
    Select positions of the pooled genome to mutate.

    Parameters
    ----------
    population : 
        instance of population class.
    L_var : int
        lenght or variant genome
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
        raise error if the lenght of active_variants is not same as 
        variants_number.

    Returns
    -------
    positions: array
        positions that will undergo substitution
    active_variants : list
        list of "active" variants, coord map onto individuals
    sub_number: int
        number of substitution events happening
    '''
    # # raise error if the vector of probabilites per site is not as long as 
    # # the variant genome
    # if L_var != len(p_sites_var):
    #     raise ValueError('p_sites_var must have L_var lenght')
    
    # list of "active" variants, coord map onto individuals
    active_variants = [[i, j] for i in population.infected_i for j 
                       in range(0,population.individuals[i]['IH_virus_number'])]
     
    # raise error if the lenght of active_variants is not same as variants_number       
    if population.active_variants_n != len(active_variants):
        raise ValueError('variants_number must have lenght "active_variants"')
    
    # vector of pooled genome positions
    positions = np.arange(0,population.active_variants_n*L_var,1)
    
    # # vector of pooled genome sites rate
    # p_sites = p_sites_var*population.active_variants_n
    # p_sites_normed = [i/np.sum(p_sites) for i in p_sites]
    
    # number of substitutions event in time step dt
    sub_number = sub_events(rng, e, dt, population.active_variants_n)
    
    if sub_number == 0:
        return 'No substitutions'
    else:    
        # return vector of positions that will be mutated p=p_sites_normed
        return rng.choice(positions,sub_number), active_variants, sub_number

def map_to_var(positions,L_var):
    '''
    Maps positions to be mutated onto active_variants.

    Parameters
    ----------
    positions: array
        positions that will undergo substitution
    L_var : int
        lenght or variant genome

    Returns
    -------
    sub_coord_varmap : tuple
        variants mutation coordinates:
            (LIST index of active_variant to mutate,
             position in the genome to mutate).
        the index here does not refer to the index of the population.individuals dic
        

    '''
    sub_coord_varmap = []
    for position in positions:
        c = [int(np.floor(position/L_var)),position%L_var]
        sub_coord_varmap.append(c)
        
    return sub_coord_varmap

def map_to_dic(active_variants, positions, L_var):
    '''
    Map mutation positions onto dictionary of population.individuals

    Parameters
    ----------
    active_variants : list
        list of "active" variants, coord map onto individuals
    positions: array
        positions that will undergo substitution 

    Returns
    -------
    sub_coord_dicmap : list
    
    coordinates of substitution:
        
    index of individual (in population.individuals) that hosts the variant    
    index of variant inside the host
    position in the genome of variant to be mutated
    
    '''
    
    sub_coord_dicmap = []
    sub_coord_varmap = map_to_var(positions, L_var)
    
    for i in sub_coord_varmap:
        # i[0] --> index of active_variants that refers to variants to mutate 
        var_coord = active_variants[i[0]]
        
        # var_coord[0] --> index of individual (in population.individuals) that hosts the variant
        # var_coord[1] --> index of variant inside the host
        # i[1] --> position in the genome of variant to be mutated
        sub_coord_dicmap.append([var_coord[0], 
                           var_coord[1],
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
        
    index of individual (in population.individuals) that hosts the variant    
    index of variant inside the host
    position in the genome of variant to be mutated
    nitrogenous base of the variant at that position

    '''
    # loop over list of coordinates to fetch the corresponding base
    for coord in sub_coord_dicmap:
        # variant --> list of mutations that a variant holds (property of individual)        
        variant = population.individuals[coord[0]]['viral_genomes'][coord[1]]
        # get mutation if the variant already have one at position
        mut_coord = [mut for mut in variant if coord[2] in mut]
        # if so, fetch the base from variant and delete mutation from variant (to add the new one)
        if mut_coord:
            coord.append(mut_coord[0][1])
            mut_index = [i for i, mut in enumerate(variant) if coord[2] in mut]
            population.individuals[coord[0]]['viral_genomes'][coord[1]].pop(mut_index[0])
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
    be included (assign the substitutions to variants with the positions in their 
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
        
    index of individual (in population.individuals) that hosts the variant    
    index of variant inside the host
    position in the genome of variant to be mutated
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

def varname(phylodots,variant_name):    
    # function to name the variants during the simulation
    n_dots = variant_name.count('.')
    try: 
        phylodots[n_dots] += 1
        variant_name = variant_name + '.' + str(phylodots[n_dots])
    except: 
        phylodots.append(0)
        variant_name = variant_name + '.0'
    
    return phylodots, variant_name

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
        parent = population.individuals[coord[0]]['IH_virus_names'][coord[1]]
        
        # add mutation to ih_lineages
        population.individuals[coord[0]]['viral_genomes'][coord[1]].append([coord[2],coord[3]])
        
        # rename the lineage after new mutation introduced
        update = varname(population.phylodots,parent)
        population.phylodots = update[0]
        child = update[1]
        population.individuals[coord[0]]['IH_virus_names'][coord[1]] = child
        
        # update lineages count
        population.individuals[coord[0]]['lineages_number'] = len(
            set(population.individuals[coord[0]]['IH_virus_names']))
        
        # save variant info for phylogenetic tree
        genome = copy.deepcopy(
            population.individuals[coord[0]]['viral_genomes'][coord[1]])
        population.track_phylogeny.append({'lineage_name': child,
                      'parent'         : parent,
                      'individual'     : coord[0],
                      'genome'         : genome,
                      'time_emergence' : population.time,
                      'host_type'      : population.individuals[coord[0]]['type'],
                      'fitness'        : dis.hamming(population.individuals[coord[0]]['viral_genomes'][coord[1]])
                    })
        
        population.lineages[child] = genome