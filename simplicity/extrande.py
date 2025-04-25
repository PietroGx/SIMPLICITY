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
from tqdm import tqdm

class ProgressReporter:
    def __init__(self, total_time, simulation_id):
        self.total_time = total_time
        self.current_time = 0
        self.step_size = 0
        
        self.leap_counter = 0
        self.reactions_counter = 0
        self.thinning_counter = 0
        self.reaction_id = None

        self.start_wall = time.time()

        self.pbar = tqdm(
                    total=total_time,
                    desc=f"RUNNING SIMPLICITY SIMULATION:  {simulation_id}",
                    unit="step",
                    bar_format="{l_bar}{bar}| {n:.2f}/{total_fmt} {unit} {postfix}",
                    position=1,
                    leave=False
                    )

    def update(self, population, delta_t, reaction_id=None, event_type=None):
        self.step_size = delta_t
        self.reaction_id = reaction_id
        self.current_time += delta_t
        self.infected = population.infected

        # Count thinning vs accepted events
        if event_type == "leap":
            self.leap_counter += 1
        elif event_type == "reaction":
            self.reactions_counter += 1
            if reaction_id == "thinning":
                self.thinning_counter += 1

        # Compute leap and thinning percentage
        total_events = self.leap_counter + self.reactions_counter
        leap_pct     =  100 * self.leap_counter / total_events if total_events > 0 else 0
        thinning_pct = 100 * self.thinning_counter / self.reactions_counter if self.reactions_counter > 0 else 0

        # Elapsed CPU time
        cpu_time = time.time() - self.start_wall

        # Update tqdm
        self.pbar.update(delta_t)
        self.pbar.set_postfix({
            "time": f"{self.current_time:.2f}",
            "step": f"{self.step_size:.2e}",
            "leap%": f"{leap_pct:.1f}%",
            "thin%": f"{thinning_pct:.1f}%",
            "last react": self.reaction_id if self.reaction_id is not None else "-",
            "infected": self.infected,
            "CPU(s)": f"{cpu_time:.1f}",
            "Time left": f"{(cpu_time / self.current_time * (self.total_time - self.current_time)):.0f}s"
        })

    def close(self):
        self.pbar.close()


def get_helpers(phenotype_model, parameters, rng1, rng2):
    """
    Returns helper functions for extrande (SIMPLICITY engine).
    """
    # model parameters
    R = parameters["R"]
    ih_tau_2 = 3.91
    tau_inf = parameters['tau_3'] + ih_tau_2
    beta = R / tau_inf
    print(f'R: {R} Beta: {beta}')
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
    
    def fire_reaction(population, propensities, tau_2):
        a0 = sum(rate for _, rate, _ in propensities)
        if tau_2 > a0:
            return 'thinning' # thinning
        threshold = 0
        for reaction_id, rate, action in propensities:
            threshold += rate
            if tau_2 <= threshold:
                action()
                return reaction_id
    
    def reaction_step(population, B):
        propensities, _ = SIDR.SIDR_propensities(population, beta, k_d, k_v, seq_rate)
        tau_2 = rng2.uniform(0, B)
        return fire_reaction(population, propensities, tau_2)
    
    def check_stop_conditions(population, t, start_time):
        if population.infected == 0:
            return "No infected left - ending simulation"
        if population.susceptibles == 0:
            return "No susceptibles left - ending simulation"
        if t > 60.0 and population.infectious_normal == 0:
            return "No infectious left - ending simulation"
        if len(population.reservoir_i) < 1000:
            return "Reservoir depleted - ending simulation"
        if time.time() - start_time > max_runtime:
            return "Max runtime exceeded"
        return None

    return {
    "compute_upperbound": compute_upperbound,
    "look_ahead": look_ahead,
    "draw_delta_t": draw_delta_t,
    "check_stop_conditions": check_stop_conditions,
    "update_step": update_step,
    "update_fitness_step": update_fitness_step,
    "mutation_step": mutation_step,
    "reaction_step": reaction_step,
    }

def extrande_core_loop(parameters, population, helpers, sim_id):
    """
    Core extrande loop.
    """
    print(f'{sim_id}')
    t = parameters['t_0']
    final_time = parameters['final_time']
    start_time = time.time()
    t_day = 0
    
    reporter = ProgressReporter(total_time=final_time, simulation_id=sim_id)
    
    min_update_threshold = 0.04  # minimum dt for intra-host update (1h step)
    dt_accumulated = 0  # initialize accumulator
    
    while t < final_time:

        L = helpers["look_ahead"](t, final_time)
        B = helpers["compute_upperbound"](population)
        
        delta_t = helpers["draw_delta_t"](B)
        dt_accumulated += delta_t
        
        reaction_id = None
        
        if delta_t > L:
            event_type = 'leap'
            # update time
            delta_t = L
            t += delta_t
            t = round(t, 10)
            population.update_time(t)
            # update system  (IH host model states)
            helpers["update_step"](population, delta_t)
            
            population.DEBUG_update_ih.append([t,L,True])
            # mutations
            helpers["mutation_step"](population, delta_t)
            # Reset dt_accumulated 
            dt_accumulated = 0
        else:
            event_type = 'reaction'
            # update time
            t += delta_t
            t = round(t, 10)
            population.update_time(t)
            # update system  (IH host model states)
            # Only update the intra-host state if accumulated dt exceeds threshold:
            if dt_accumulated >= min_update_threshold:
                helpers["update_step"](population, dt_accumulated)
                population.DEBUG_update_ih.append([t,dt_accumulated,False])
                # mutations
                helpers["mutation_step"](population, dt_accumulated)
                # update fitness
                helpers["update_fitness_step"](population)
                dt_accumulated = 0  # reset accumulator
            # reactions
            reaction_id = helpers["reaction_step"](population, B)
        # update progress bar
        reporter.update(population, delta_t, reaction_id, event_type)

        # update system trajectory
        population.update_trajectory()
        if math.floor(t) > t_day:
            t_day += 1
            population.update_lineage_frequency_t(t)
        # break if conditions met
        reason = helpers["check_stop_conditions"](population, t, start_time)
        if reason:
            break

    reporter.close()
    if reason:
        tqdm.write(f"\n {reason}")
    
    tqdm.write("\n----------------------------------------")
    tqdm.write("SIMPLICITY SIMULATION COMPLETED")
    tqdm.write("----------------------------------------\n")
    tqdm.write('')
    return population

def extrande_factory(phenotype_model, parameters, sim_id, rng1, rng2):
    
    def extrande_generic(population):
        helpers = get_helpers(phenotype_model, parameters, rng1, rng2)
        return extrande_core_loop(parameters, population, helpers, sim_id)

    return extrande_generic
