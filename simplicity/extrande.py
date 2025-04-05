#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 28 13:31:54 2025

@author: pietro
"""
import time
import math
import numpy as np
import simplicity.population_model as SIDR
import simplicity.evolution.mutations as evo
import simplicity.evolution.reference as ref
import simplicity.phenotype.consensus as c
import simplicity.phenotype.update as pheno
import simplicity.tuning.diagnosis_rate as dr

def get_helpers(phenotype_model, parameters, rng1, rng2):
    """
    Returns helper functions for extrande (SIMPLICITY engine).
    """
    # model parameters
    R = parameters["R"]
    ih_tau_2 = 3.91
    tau_inf = parameters['tau_3'] + ih_tau_2
    beta = R / tau_inf
    k_d = dr.get_k_d_from_diagnosis_rate(parameters["diagnosis_rate"], parameters["tau_3"])
    k_v = parameters["IH_virus_emergence_rate"]
    e = np.sum([parameters["nucleotide_substitution_rate"]] * len(ref.get_reference()))
    seq_rate = parameters["sequencing_rate"]
    max_runtime = parameters["max_runtime"]
    # for fitness update
    update_all_fitness = pheno.update_fitness_factory(phenotype_model)
    use_consensus = phenotype_model == "immune_waning"
    last_consensus_snapshot = {"t_snapshot": 0, "consensus": []} if use_consensus else None
    
    def compute_upperbound(population):
        propensities, params = SIDR.SIDR_propensities(population, beta, k_d, k_v, seq_rate)
        # a0 = sum(rate for rate, _ in propensities)
        B = np.sum(params) * population.infected
        # print(f"B = {B}")
        return B
    
    def look_ahead(t, final_time):
        '''
        We clamp the look ahead to 1 day: if extrande does a step bigger than a day
        we need to update the ih model to make sure not to miss host transitions 
        and skewer their intra-host trajectory
        '''
        # L =  min(final_time - t, population.get_next_ih_transition() - t)
        # return max(L, 0.0417) # either next IH state update or 1/24 = 1 hour
        return min(1, final_time - t) 

    def draw_delta_t(B):
        tau = rng1.uniform(0, 1)
        delta_t = -np.log(tau) / B
        delta_t = round(delta_t, 10)
        return delta_t

    def update_step(population, delta_t):
        population.update_states(delta_t)
    
    def update_fitness_step(population):
        t = population.time
        individuals_to_update = sorted(population.infected_i)
        # calculate consensus every 5 simulation days and use it for fitness score (immune waning model)
        if use_consensus:
            if np.floor(t) > last_consensus_snapshot["t_snapshot"]:
                last_consensus_snapshot["t_snapshot"] += 5
                consensus = c.get_consensus(population.consensus_snapshot, t)
                population.consensus_sequences_t.append([consensus, t])
                last_consensus_snapshot['consensus'] = consensus 
            else:
                consensus = last_consensus_snapshot['consensus']
            update_all_fitness(population, individuals_to_update, consensus)
        # fitness update for linear model
        else:
            update_all_fitness(population, individuals_to_update)
        
        population.update_fitness_trajectory()
    
    def mutation_step(population, delta_t):
        delta_t_y = delta_t / 365.25 # time in years
        if use_consensus:
            consensus = last_consensus_snapshot['consensus']
            evo.mutate(population, e, delta_t_y, phenotype_model, consensus)
        else:
            evo.mutate(population, e, delta_t_y, phenotype_model)
    
    def fire_reaction(population, reactions, tau_2):
        a0 = sum(rate for rate, _ in reactions)
        if tau_2 > a0:
            return  1 # thinning
        threshold = 0
        for rate, action in reactions:
            threshold += rate
            if tau_2 <= threshold:
                action()
                return 0
    
    def reaction_step(population, B):
        reactions, _ = SIDR.SIDR_propensities(population, beta, k_d, k_v, seq_rate)
        tau_2 = rng2.uniform(0, B)
        return fire_reaction(population, reactions, tau_2)

    def report_progress(t, final_time, last_progress):
        progress = t / final_time
        if progress >= last_progress + 0.1:
            print(f"Progress: {progress*100:.0f}% (t = {t:.2f})")
            return progress
        return last_progress

    def check_stop_conditions(population, t, start_time):
        
        if population.infected == 0:
            print("\n No infected left - ending simulation")
            return True
        if population.susceptibles == 0:
            print("\nNo susceptibles left - ending simulation")
            return True
        if t > 60.0:
            if population.infectious_normal == 0:
                print("\n No infectious left - ending simulation")
                return True
        if len(population.reservoir_i) < 1000:
            print("\n Reservoir depleted - ending simulation")
            return True
        if time.time() - start_time > max_runtime:
            print("\n Max runtime exceeded")
            return True
        return False

    return {
    "report_progress": report_progress,
    "compute_upperbound": compute_upperbound,
    "look_ahead": look_ahead,
    "draw_delta_t": draw_delta_t,
    "check_stop_conditions": check_stop_conditions,
    "update_step": update_step,
    "update_fitness_step": update_fitness_step,
    "mutation_step": mutation_step,
    "reaction_step": reaction_step,
    }

def extrande_core_loop(parameters, population, helpers):
    """
    Core extrande loop.
    """
    t = parameters['t_0']
    final_time = parameters['final_time']
    start_time = time.time()
    t_day = 0
    last_progress = 0
    
    leap_counter = 0
    extrande_counter = 0
    thinning_counter = 0
    
    min_update_threshold = 0.04  # minimum dt for intra-host update (1h step)
    dt_accumulated = 0  # initialize accumulator
    
    while t < final_time:
        last_progress = helpers["report_progress"](t, final_time, last_progress)

        L = helpers["look_ahead"](t, final_time)
        B = helpers["compute_upperbound"](population)
        # print(f'Upper bound: {B}')
        delta_t = helpers["draw_delta_t"](B)
        dt_accumulated += delta_t
        
        if delta_t > L:
            leap_counter +=1
            # update time
            t += L
            t = round(t, 10)
            population.update_time(t)
            # update system  (IH host model states)
            helpers["update_step"](population, L)
            # mutations
            helpers["mutation_step"](population, L)
            # Reset dt_accumulated 
            dt_accumulated = 0
        else:
            extrande_counter+=1
            # update time
            t += delta_t
            t = round(t, 10)
            population.update_time(t)
            # update system  (IH host model states)
            # Only update the intra-host state if accumulated dt exceeds threshold:
            if dt_accumulated >= min_update_threshold:
                helpers["update_step"](population, dt_accumulated)
                # mutations
                helpers["mutation_step"](population, dt_accumulated)
                # update fitness
                helpers["update_fitness_step"](population)
                dt_accumulated = 0  # reset accumulator
            # reactions
            thinning_counter += helpers["reaction_step"](population, B)
            
        # print(t)
        # print(f"L = {L}, delta_t = {delta_t}")

        # update system trajectory
        population.update_trajectory()
        if math.floor(t) > t_day:
            t_day += 1
            population.update_lineage_frequency_t(t)
        # break if conditions met
        if helpers["check_stop_conditions"](population, t, start_time):
            break

    print("\n----------------------------------------")
    print("SIMPLICITY SIMULATION COMPLETED")
    print("----------------------------------------\n")
    # print(f'Extrande counter: {extrande_counter}')
    # print(f'Leap counter: {leap_counter}')
    print(f'Leaps Ratio (leaps/total steps): {leap_counter/(extrande_counter+leap_counter)}')
    print(f'Thinning Ratio (thinning/total reactions firing): {thinning_counter/extrande_counter}')
    print("----------------------------------------\n")
    return population

def extrande_factory(phenotype_model, parameters, rng1, rng2):
    
    def extrande_generic(population):
        helpers = get_helpers(phenotype_model, parameters, rng1, rng2)
        return extrande_core_loop(parameters, population, helpers)

    return extrande_generic
