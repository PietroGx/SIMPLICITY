"""
EXTRANDE implementation of the SIR evolutionary model of
SARS-CoV-2 infection dynamics

@author: Pietro Gerletti
"""
import time 
import math
import numpy as np
import simplicity.phenotype.consensus as c
import simplicity.phenotype.update as pheno
import simplicity.evolution.reference as ref
import simplicity.tuning.diagnosis_rate as dr
# from memory_profiler import profile
    
def extrande_factory(phenotype_model, parameters,
                     rng1, rng2, path):
    '''
    Returns
    -------
    selected versin of extrande function

    '''
    # check if label was passed correctly
    if phenotype_model == 'distance from wt' or phenotype_model == 'immune waning':
        pass
    else:
        raise ValueError('check the arguments of the factory')
    
    def beta(R,avg_infectious_time):
        # return infection rate
        return R / avg_infectious_time
    
    def upperbound(beta_normal,k_d,k_v, population):
        # compute upper bound for extrande algorithm 
        a1 = beta_normal*population.infectious_normal
        a2 = k_d*population.detectables 
        a3 = k_v*population.infected
        a0 = a1 + a2 + a3 
        return a0*parameters["F"]

    # time limits
    final_time = parameters['final_time']
    max_runtime = parameters['max_runtime']
    
    # starting value for look-ahead time horizon
    L0 = 1
    # virus reproduction number for individuals in the population
    R = parameters["R"]
    # average infectious time for individuals infected with SARS-CoV-2
    ih_tau_2 = 3.91
    tau_inf = parameters['tau_3'] + ih_tau_2
    # diagnosis rate
    k_d = dr.get_k_d_from_diagnosis_rate(parameters["diagnosis_rate"])
    # multiple variants rate
    k_v = parameters["IH_virus_emergence_rate"]
    # per site per year SPIKE substitution rate
    avg_spike_ssy = parameters["evolutionary_rate"]
    # substitution rates (per site, pe_r year) of SARS-CoV-2 Spike gene
    rate_sites_var = [avg_spike_ssy for i in range(len(ref.get_reference()))]
    # substitution rate (per year) of SARS-CoV-2 variants
    e = np.sum(rate_sites_var)
    # set sequencing rate
    seq_rate = parameters["sequencing_rate"]

    # extrande with simple phenotype model
    # @profile
    def extrande_0(population):
        # count number of infectious and detectables individuals
        population.update_inf_det()
        # record computation start time
        start00 = time.time()
        # start time
        t = parameters['t_0']
        
        # count time leaps and thinning reactions
        # thinning = 0
        # leap = 0
        # compute beta
        b = beta(R,tau_inf)
        
        # counter to update lineages count
        t_day = 0 # only count lineages frequency once per day
        
        # store times of infections and diagnosis reaction (benchmarking)
        # store_inf_t  = []
        # store_dia_t  = []
        last_printed_progress = 0 
        ###########################################################################
        # repeat (EXTRANDE loop)
        while t<final_time:
            
            # Calculate progress as a fraction of t_final
            progress = t / final_time
            if progress >= last_printed_progress + 0.1:  # 10% increments
                print(f"Progress: {progress*100:.0f}% (t/final_time = {t:.2f}/{final_time} d)")
                last_printed_progress = progress  # Update the last printed progress
            
            # look ahead time horizon
            L = L0/population.infected
            
            # compute upper bound 
            B = upperbound(b,k_d,k_v, population)
           
            ## select timepoint of next reaction
            # draw tau
            tau_1 = rng1.uniform(0,1)
            
            # compute delta_t from tau_1 (equivalent of exp distributed delta_t)
            delta_t = -1/B * np.log(tau_1)
            
            
                
            # time leap
            #######################################################################
            if delta_t > L: 
                # update time 
                t += L
                # count leaps
                # leap += 1
            
            # reaction happens 
            #######################################################################
            else:
                # update time 
                t += delta_t
                
                ###############################################################
                population.update_time(t)
                # update states 
                population.update_states(t-population.trajectory[-1][0])
                
                # remove subjects that are not infectious anymore 
                population.recovery()
                
                # update the number of detectable and infectious individuals
                population.update_inf_det()
                ###############################################################
                
                # mutations happen
                delta_t_y = delta_t/365.25 # delta t in years
                population.mutate(e, delta_t_y, phenotype_model,'')
                
                # reaction 1 - infection
                a1 = b*population.infectious_normal
                # reaction 2 - detection
                a2 = k_d*population.detectables
                # reaction 3 - multiple variants
                a3 = k_v*population.infected
                # compute a0
                a0 = a1+a2+a3
               
                # reactions happen
                tau_2 = rng2.uniform(0,B)
                         
                if tau_2 > a0:     # thinning reaction
                    # thinning += 1
                    pass
                else:
                    if a1+a2 < tau_2  <= a0:  # reaction 3 - virus duplication
                        population.add_variant()
                    elif a1 < tau_2  <= a1+a2:  # reaction 2 - diagnosis
                        population.diagnosis(seq_rate)
                        
                        # store_dia_t.append(t)
                                            
                    elif tau_2 <= a1:           # reaction 1 - infection
                        population.infection()
                        
                        # store_inf_t.append(t)
                        
                    # save the system state at time t
                    population.update_trajectory()
                    # update count of lineages frequency
                    if math.floor(t) > t_day: # only count once per day
                        t_day += 1    
                        population.count_variants_t(t)
                
                ##  break if conditions are met
                if population.infected == 0:
                    print('')
                    print('no infected left')
                    break
                if population.susceptibles == 0:
                    print('')
                    print('No susceptibles left - ending simulation')
                    break
                if population.infectious_normal == 0:
                    print('')
                    print('no infectious left')
                    break
                if len(population.reservoir_i) < 1000:
                    print('')
                    print('reservoir depleted - ending simulation')
                    break
                if time.time()-start00 > max_runtime:
                    print('')
                    print('Elapsed time exceeds limit')
                    break
        print('')     
        print('----------------------------------------')
        # print('thinning: ', thinning)
        # print('Total reactions: ', len(population.trajectory))
        # print('Thinning/total reactions',thinning/(len(population.trajectory)))
        print('SIMPLICITY SIMULATION COMPLETED')
        print('----------------------------------------')
        print('')
        # Save lists to text files 
        # with open(path+'/diagnosis_times.csv', 'w') as f:
        #     for item in store_dia_t:
        #         f.write(f"{item}\n")
        # with open(path+'/infection_times.csv', 'w') as f:
        #     for item in store_inf_t:
        #         f.write(f"{item}\n")
        
        return population # population after simulation is called simulation_output
    
    # extrande with immune-waning phenotype model
    # @profile
    def extrande_1(population):
        # count number of infectious and detectables individuals
        population.update_inf_det()
        
        # record computation start time
        start00 = time.time()
        # start time
        t = parameters['t_0']
       
        # count time leaps and thinning reactions
        # thinning = 0
        # leap = 0
        # compute beta
        beta_normal = beta(R,tau_inf)
        # counter to update lineages count
        t_day = 0 # only count lineages frequency once per day
        
        t_snapshot = 0
        consensus  = []
        
        # store times of infections and diagnosis reaction (benchmarking)
        # store_inf_t  = []
        # store_dia_t  = []
        # pop_size_runtime_plot_coord = []
        last_printed_progress = 0 
        ###########################################################################
        # repeat (EXTRANDE loop)
        while t<final_time:
            
            # Calculate progress as a fraction of t_final
            progress = t / final_time
            if progress >= last_printed_progress + 0.1:  # 10% increments
                print(f"Progress: {progress*100:.0f}% (t/final_time = {t:.2f}/{final_time} d)")
                last_printed_progress = progress  # Update the last printed progress
           
            L  = L0/population.infected
           
            # compute upper bound 
            B = upperbound(beta_normal,k_d,k_v, population)
          
            ## select timepoint of next reaction
            # draw tau
            tau_1 = rng1.uniform(0,1)
           
            # compute delta_t from tau_1 (equivalent of exp distributed delta_t)
            delta_t = -1/B * np.log(tau_1)
            
            # time leap
            #######################################################################
            if delta_t > L: 
                # update time 
                t += L
                # count leaps
                # leap += 1                
            
            # reaction happens 
            #######################################################################
            else:
                # update time 
                t += delta_t
                
                ###############################################################
                population.update_time(t)
                # update states 
                population.update_states(t-population.trajectory[-1][0])
                
                # remove subjects that are not infectious anymore 
                population.recovery()
                
                # update the number of detectable and infectious individuals
                population.update_inf_det()
                ###############################################################
                
                # calculate consensus sequence every 10 days
                if np.floor(t) > t_snapshot:
                    t_snapshot += 5
                    consensus =c.get_consensus(population.consensus_snapshot,t)
                        
                # update fitness values of all individuals in the simulation
                individuals_to_update = population.infected_i
                update_all_fitness = pheno.update_fitness_factory(phenotype_model)
                population.individuals = update_all_fitness(
                    population.individuals, individuals_to_update, consensus)
                
                # store mean fitness for plot
                population.track_fitness()
                
                # mutations happen
                delta_t_y = delta_t/365.25 # delta t in years
                population.mutate(e, delta_t_y, phenotype_model, consensus)
                
                # reaction 1 - infection
                a1 = beta_normal*population.infectious_normal
                # reaction 2 - detection
                a2 = k_d*population.detectables
                # reaction 3 - multiple variants
                a3 = k_v*population.infected
                # compute a0
                a0 = a1+a2+a3
                
                # reactions happen
                tau_2 = rng2.uniform(0,B)
                         
                if tau_2 > a0:     # thinning reaction
                    # thinning += 1
                    pass
                
                else:
                    if a1+a2 < tau_2  <= a0:  # reaction 3 - virus duplication
                        population.add_variant()
                    elif a1 < tau_2  <= a1+a2:  # reaction 2 - diagnosis
                        population.diagnosis(seq_rate)
                        # store_dia_t.append(t)
                    elif tau_2 <= a1:           # reaction 1 - infection
                        population.infection()
                        # store_inf_t.append(t)
                    # save the system state at time t
                    population.update_trajectory()
                
                
                # pop_size_runtime_plot_coord.append((time.time()-start00, population.infected))
                
                # update count of lineages frequency
                if math.floor(t) > t_day: # only count once per day
                    t_day += 1    
                    population.count_variants_t(t)
            
                ##  break if conditions are met
                if population.infected == 0:
                    print('')
                    print('no infected left - ending simulation')
                    break
                if population.susceptibles == 0:
                    print('No susceptibles left - ending simulation')
                    break
                if population.infectious_normal == 0:
                    print('')
                    print('no infectious left - ending simulation')
                    break
                if len(population.reservoir_i) < 1000:
                    print('')
                    print('reservoir depleted - ending simulation')
                    break
                if time.time()-start00 > max_runtime:
                    t = final_time
                    print('')
                    print('Elapsed time exceeds limit')
                    break
        print('')     
        print('----------------------------------------')
        # print('thinning: ', thinning)
        # print('Total reactions: ', len(population.trajectory))
        # print('Thinning/total reactions',thinning/(len(population.trajectory)))
        print('SIMPLICITY SIMULATION COMPLETED')
        print('----------------------------------------')
        print('')
        # # Save lists to text files 
        # with open(path+'/diagnosis_times.csv', 'w') as f:
        #     for item in store_dia_t:
        #         f.write(f"{item}\n")
        
        # with open(path+'/infection_times.csv', 'w') as f:
        #     for item in store_inf_t:
        #         f.write(f"{item}\n")
       
        # import csv
        # with open('extrande_pop_runtime.csv', mode='w', newline='') as file:
        #     writer = csv.writer(file)
        #     writer.writerows(pop_size_runtime_plot_coord)
            
        return population # population after simulation is called simulation_output
    
    if phenotype_model == 'distance from wt':
        # phenotype model: distance from wt
        # return extrande with simple phenotype model
        return extrande_0
        
    elif phenotype_model == 'immune waning':
        # phenotype model: immune waning    
        # return extrande with immune-waning phenotype model
        return extrande_1
   

    
    
       