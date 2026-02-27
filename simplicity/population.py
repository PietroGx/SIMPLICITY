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

"""
Created on Tue Jun  6 13:13:14 2023

@author: pietro
"""
import simplicity.intra_host_model    as h
import simplicity.evolution.reference as ref
from   simplicity.random_gen          import randomgen
import pandas as pd
import numpy as np
import scipy.stats

class Population:
    '''
    The class defines a population for the SIMPLICITY simulations. It contains
    the data about every individual as well as their intra-host model.
    '''
    def __init__(self,
                 size,I_0,
                 ih_model_parameters,
                 rng3,rng4,rng5,rng6,
                 NSR_long,
                 long_shedders_ratio=0,
                 sequence_long_shedders=False,
                 reservoir=100000):
        
        # random number generator
        self.rng3 = rng3 # for intra-host model states update
        self.rng4 = rng4 # for electing individuals|lineages when reactions happen
        self.rng5 = rng5 # for mutation model
        self.rng6 = rng6 # for synthetic sequencing data
        
        self.size = size
        self.long_shedders_ratio = long_shedders_ratio
        
        # compartments and population groups ---------------------------------
        
        self.susceptibles = size - I_0 # number of susceptible individuals
        self.infected     = I_0        # number of infected individuals
        self.diagnosed    = 0          # number of diagnosed individuals
        self.recovered    = 0          # number of recovered individuals

        self.long_shedders = 0   # number of long shedders
        
        self.infectious_standard = 0
        self.infectious_long     = 0
        self.infectious          = 0 # number of infectious individuals
        
        self.detectables_standard = 0
        self.detectables_long     = 0
        self.detectables          = 0 # number of detectable individuals
        
        self.reservoir = reservoir   # size of total population (not everyone is 
                                     # susceptible at the beginning, when 
                                     # individuals get removed from the system
                                     # new ones from the reservoir become 
                                     # susceptible) 
                                     
        # individuals ---------------------------------------------------------
        self.individuals = {}          # store individuals data
        self.reservoir_i    = set()    # set of indices of individuals in the reservoir
        
        self.susceptibles_i = set()    # set of susceptible individuals indices  
        self.infected_i     = set()    # set of infected individuals indices 
        self.diagnosed_i    = set()    # set of diagnosed individuals indices
        self.recovered_i    = set()    # set of recovered individuals indices
        
        self.long_shedder_i = set()    # set of long shedder individuals
        
        
        self.infectious_standard_i = set()
        self.infectious_long_i     = set()
        self.infectious_i          = set() # set of infectious individuals indices
        
        self.detectables_standard_i = set()
        self.detectables_long_i     = set()
        self.detectables_i         = set() # set of detectable individuals indices 
        
        self.exclude_i      = set()    # set to store the newly infected individual (excludes it from states update at time of infection)
        
        # set attributes for evolutionary model ===============================
        self.ref_genome      = ref.get_reference()  # sequence of reference genome
        self.L_lin           = len(self.ref_genome) # lenght of reference genome
        self.active_lineages_n = I_0                # number of IH viruses in the 
                                                    # infected population
        self.NSR_long = NSR_long # long shedders substitution rate
        
        # Phylogenetic tree data ----------------------------------------------
        self.sequence_long_shedders = sequence_long_shedders
        self.phylogenetic_data = [{  'Time_emergence'  : 0,
                                     'Lineage_name'    : 'wt',
                                     'Lineage_parent'  : None,
                                     'Genome'          : [],
                                     'Host_type'       : 'standard',
                                     'Total_infections': 0
                                 }]
        self._phylo_name_map = { row['Lineage_name']: row for row in self.phylogenetic_data }
        self.phylodots = []         # needed to name lineages
        self.sequencing_data = []   # to store sequencing data
        # ---------------------------------------------------------------------
        
        self.lineage_frequency = [] # count lineage frequency in the population
        
        # -------------------------------------------------------------------------   

        self.consensus_snapshot = [[[],1,0]] # [sequence, lineage_frequency, w_t(t_sim)]
        self.consensus_sequences_t = [[[],0]] # list to store consensus everytime is calculated in a simulation [consensus,t]
        
        # system trajectory ---------------------------------------------------
        self.time         = 0 
        self.trajectory = [[self.time,
                           self.infected,
                           self.diagnosed,
                           self.recovered,
                           self.infectious,
                           self.detectables,
                           self.susceptibles,
                           self.long_shedders
                           ]]
        
        self.last_infection = {}    # tracks the information about the last infection event
                                    # that happened in the simulaiton, used for 
                                    # R_effective calculations
        self.R_effective_trajectory = []
        self.fitness_trajectory = []
        
       # ih model -------------------------------------------------------------
        self.update_ih_mode = 'matrix'
        self.host_model = {'standard': 
                                  h.Host(tau_1=ih_model_parameters["tau_1"],
                                  tau_2=ih_model_parameters["tau_2"],
                                  tau_3=ih_model_parameters["tau_3"],
                                  tau_4=ih_model_parameters["tau_4"],
                                  update_mode = self.update_ih_mode) ,  # intra host model for standard individuals 
        
                           'long_shedder': 
                                  h.Host(tau_1=ih_model_parameters["tau_1"],
                                  tau_2=ih_model_parameters["tau_2"],
                                  tau_3=ih_model_parameters["tau_3_long"],
                                  tau_4=ih_model_parameters["tau_4"],
                                  update_mode = self.update_ih_mode)   # intra host model for long-shedders
                           }
        # -------------------------------------------------------------------------
        # counter for inf reactions
        self.count_infections = 0
        # self.count_infections_from_long_shedders = 0
        
        # dictionary with all individuals data 
        self.individuals = self._init_individuals(size,I_0)
        
        
    # -------------------------------------------------------------------------   
    # -------------------------------------------------------------------------
    def _init_individuals(self,size,I_0):
        '''
        Create dictionary with all individuals in the simulation

        Parameters
        ----------
        size : int
            Size of the population.
        I_0 : int
            Number of infected individuals at the beginning of the simulation.

        Returns
        -------
        dic : dict
            Dictionary containing the info of all individuals in the population.

        '''
        
        dic = {}
        # create an entry in the dictionary for each individual (0 to number of
        # total individuals in the population (reservoir))
        for i in range(self.reservoir):
            
            dic[i] = {
                     't_infection' : None,
                     't_not_infected': None,
                     
                     't_infectious': None,
                     't_not_infectious': None,
                     
                     'type'        : 'standard',
                     'state_t'     : 0,
                     't_next_state': None,
                     'state'       : 'susceptible',
                     
                     'parent'      : None,
                     'inherited_lineage': None,
                     'new_infections'  : [],
                     
                     'IH_lineages'   : [],
                     'IH_unique_lineages_number': 0, #1
                     'IH_lineages_number'    : 0,    #1
                     'IH_lineages_max': self.rng3.integers(1,5),
                     'IH_lineages_fitness_score' : [], #1
                     'IH_lineages_trajectory': {}, # lineage name : [ih_time_start, ih_time_end]
                     'time_last_weight_event': 0, # time since infection or last mutation
                     'fitness_score'     :  1e-6  # fitness floor (individual fitness)
                    }
            # add index of individuals to either susceptibles indices or to 
            # reservoir indices (the simulaiton starts with a pool of susceptibles)
            # that is replenished from the reservoir every time individuals
            # are removed from the system
            if i < size:
                self.susceptibles_i.add(i) 
            else:
                self.reservoir_i.add(i)
            
        # set individuals infected at the beginning of the simulation
        for i in range(I_0): # update data of individuals infected at the beginning of the simulation
            
            dic[i]['parent']       = 'root'
            dic[i]['t_infection']  = 0
            dic[i]['t_infectious'] = None
            dic[i]['state_t']      = 0
            # sample next jump time from exp dist.
            state_t = dic[i]['state_t']
            rate = - self.host_model['standard'].A[state_t][state_t]
            dic[i]['t_next_state'] = self.rng3.exponential(scale=1/rate) 
            
            dic[i]['state']                     = 'infected'
            dic[i]['IH_lineages']               = ['wt']
            dic[i]['IH_lineages_fitness_score'] = [1e-6]
            dic[i]['IH_unique_lineages_number'] = 1 
            dic[i]['IH_lineages_number']        = 1
            dic[i]['inherited_lineage']         = 'wt'
            dic[i]['IH_lineages_trajectory']['wt'] = {'ih_birth':None,'ih_death':None}
            
            self.susceptibles_i.remove(i)  
            self.infected_i.add(i)

        # update lineage_frequency
        self.lineage_frequency.append({'Lineage_name'              :'wt',
                                       'Time_sampling'             :0,
                                       'Frequency_at_t'            :1,
                                       'Individuals_infected_at_t' : I_0
                                       })
        
        # return dictionary containing all individuals data (self.individuals)
        return dic        
    
    # -------------------------------------------------------------------------
    def get_lineage_genome(self, lineage_name):
       '''
       Fetch lineage genome from lineage name
       '''
       genome = next((d['Genome'] for d in self.phylogenetic_data if d['Lineage_name'] == lineage_name))
       # print(f'{lineage_name} Genome: ', genome)
       return genome
    # -------------------------------------------------------------------------
    #                               Updates
    # -------------------------------------------------------------------------
    
    def update_time(self,time):
        # update the time 
        self.time = time
        
    def _update_states_matrix(self, delta_t, individual_type):
        """
        General intra-host state updater for a pop group using its host_model.
        """
        
        if individual_type == 'standard':
            
            infected_to_update = [i for i in self.infected_i if i not in self.exclude_i and 
                                                                i not in self.long_shedder_i]
            
        elif individual_type == 'long_shedder':
            
            infected_to_update = [i for i in self.infected_i if i not in self.exclude_i and 
                                                                i in self.long_shedder_i]
        
        else:
            raise ValueError('Invalid individual type!')
            
        # compute transition probabilitiy vectors
        host_model = self.host_model[individual_type]
        all_probabilities = host_model.compute_all_probabilities(delta_t)
        
        # draw random variables for each infected individual in the subpopulation
        taus = self.rng3.uniform(size=len(infected_to_update))
    
        for idx, i in enumerate(infected_to_update):
            ind = self.individuals[i]
            state = ind['state_t']
            prob = all_probabilities[state]
            new_state = host_model.update_state(prob, taus[idx])
            ind['state_t'] = new_state
    
            # Recovery check
            if new_state == 20:
                ind['state'] = 'recovered'
                if ind['t_not_infectious'] is None:
                    ind['t_not_infectious'] = self.time
                if ind['t_not_infected'] is None:
                    ind['t_not_infected'] = self.time
                else:
                    raise ValueError('Individual already recovered!!')
    
                self.infected_i.remove(i)
                self.infectious_i.discard(i)
                self.detectables_i.discard(i)
                self.recovered_i.add(i)
                
                self.long_shedder_i.discard(i)
    
                self.infected -= 1
                self.recovered += 1
                self.susceptibles += 1
    
                self.active_lineages_n -= ind['IH_lineages_number']
    
                new_susceptible = self.reservoir_i.pop()
                self.susceptibles_i.add(new_susceptible)
                continue
    
            elif new_state <= 4:
                continue
    
            # Detectable update (5–19)
            if 4 < new_state < 20:
                self.detectables_i.add(i)
            else:
                self.detectables_i.discard(i)
    
            # Infectious update (5–18 for standard individuals)
            if 4 < new_state < 19:
                if i not in self.infectious_i:
                    self.infectious_i.add(i)
                    if ind['t_infectious'] is None:
                        ind['t_infectious'] = self.time
                    else:
                        raise ValueError(f'Individual {i} t_infectious already set!!')
            else:
                if new_state != 19:
                    raise ValueError('State here should be 19 only')
                if i in self.infectious_i:
                    self.infectious_i.remove(i)
                    if ind['t_not_infectious'] is None:
                        ind['t_not_infectious'] = self.time
                    else:
                        raise ValueError(f'Individual {i} t_not_infectious already set!!')
    
        # Post-processing
        self.exclude_i = set()
        
        # Update compartments
        self.infectious_standard_i = sorted(self.infectious_i - self.long_shedder_i)
        self.infectious_long_i = sorted(self.infectious_i & self.long_shedder_i)
        
        self.detectables_standard_i = sorted(self.detectables_i - self.long_shedder_i)
        self.detectables_long_i = sorted(self.detectables_i & self.long_shedder_i)
        
        self.infectious = len(self.infectious_i)
        self.infectious_standard = len(self.infectious_standard_i)
        self.infectious_long  = len(self.infectious_long_i)
        
        self.detectables = len(self.detectables_i)
        self.detectables_long = len(self.detectables_long_i)
        self.detectables_standard = len(self.detectables_standard_i)
        
        self.long_shedders = len(self.long_shedder_i)

    def update_states(self, delta_t):
        if self.update_ih_mode == "jump":
            self._update_states_jump()
        elif self.update_ih_mode == "matrix":
            self._update_states_matrix(delta_t,'standard')
            self._update_states_matrix(delta_t, 'long_shedder')
        else:
            raise ValueError(f"Unknown update_mode: {self.update_mode}")
    
    def update_trajectory(self):
        # update the system trajectory
        self.trajectory.append([self.time,
                           self.infected,
                           self.diagnosed,
                           self.recovered,
                           self.infectious,
                           self.detectables,
                           self.susceptibles,
                           self.long_shedders
                           ])
    
    def update_fitness_trajectory(self):
        fitness_scores = [self.individuals[i]['fitness_score'] for i in self.infected_i]
        if not fitness_scores:
            self.fitness_trajectory.append([self.time, 0, 0, 0])
            return
    
        mean_fitness = np.mean(fitness_scores)
        std_fitness = np.std(fitness_scores)
    
        # Entropy of normalized fitness 
        total = sum(fitness_scores)
        normed = [f / total for f in fitness_scores]
        entropy = scipy.stats.entropy(normed) if total > 0 else 0.0
    
        self.fitness_trajectory.append({
                                        'Time': self.time,
                                        'Mean': mean_fitness,
                                        'Std': std_fitness,
                                        'Entropy': entropy
                                    })
    
    def update_lineage_frequency_t(self, t):
        '''
        Count how many individuals are infected by each lineage at time
        t in the population and store the following information:
        time, lineage_name, frequency.
        '''
        # store lineages here: {'lineage name': number of infected individuals}
        count_lineages_t = {}
        
        # Loop through the infected individuals and get a list of unique IH_virus names (lineage)
        for individual_index in set(self.infected_i) - set(self.long_shedder_i):
            unique_lineages = set(self.individuals[individual_index]['IH_lineages'])
            for lineage_name in unique_lineages:
                count_lineages_t[lineage_name] = count_lineages_t.get(lineage_name, 0) + 1
    
        # Calculate the total count for normalization
        infected_individuals_at_t_total = sum(count_lineages_t.values()) # note that if an individual is infected by more than one lineage it  will count double
        
        # Store: [time, lineage_name, frequency, infected_individuals_at_t]
        for lineage_name in count_lineages_t:
            # Compute frequency as the relative count
            infected_individuals_lin = count_lineages_t[lineage_name]
            frequency = infected_individuals_lin / infected_individuals_at_t_total
            # store frequency data
            dic = {'Lineage_name'              :lineage_name,
                   'Time_sampling'             :t,
                   'Frequency_at_t'            :frequency,
                   'Individuals_infected_at_t' :infected_individuals_lin
                }
            self.lineage_frequency.append(dic)
        
        for lineage_name in count_lineages_t:
            frequency = count_lineages_t[lineage_name] / infected_individuals_at_t_total
            self.consensus_snapshot.append([self.get_lineage_genome(lineage_name),
                                                frequency,
                                                t])
            
# -----------------------------------------------------------------------------

    def update_ih_lineages_trajectories(self):
        individuals_data = pd.DataFrame(self.individuals).transpose()
        for idx, row in individuals_data.iterrows():
            lineage_traj_dic = row['IH_lineages_trajectory']  
        
            for lineage in lineage_traj_dic:
                if lineage_traj_dic[lineage].get('ih_birth') is None:
                    lineage_traj_dic[lineage]['ih_birth'] = row['t_infectious']
                if lineage_traj_dic[lineage].get('ih_death') is None:
                    lineage_traj_dic[lineage]['ih_death'] = row['t_not_infectious']
        
            individuals_data.at[idx, 'IH_lineages_trajectory'] = lineage_traj_dic
        return individuals_data
        
    def individuals_data_to_df(self):
        # return population dictionary as data frame
        df = self.update_ih_lineages_trajectories()
        filtered_df = df[~df['state'].str.contains('susceptible')]
        filtered_df = filtered_df.drop('t_next_state', axis=1)
        return filtered_df
    
    def phylogenetic_data_to_df(self):
        # return phylogeny dictionary as data frame
        return pd.DataFrame(self.phylogenetic_data)
    
    def lineage_frequency_to_df(self):
        # return lineage_frequency as data frame
        return pd.DataFrame(self.lineage_frequency)
        
    def fitness_trajectory_to_df(self):
        return pd.DataFrame(self.fitness_trajectory)

    
# -----------------------------------------------------------------------------
# =============================================================================
# -----------------------------------------------------------------------------  

def create_population(parameters):
    '''
    Create population instance from parameters file and return it.
    '''
    # population parameters
    pop_size = parameters['population_size']
    I_0      = parameters['infected_individuals_at_start']
    seed     = parameters['seed']
    
    long_shedders_ratio = parameters['long_shedders_ratio']
    M_nsr_long    = parameters['M_nsr_long']
    sequence_long_shedders = parameters['sequence_long_shedders']
    
    NSR_long = parameters['nucleotide_substitution_rate'] * M_nsr_long
    
    ih_model_parameters = {
        'tau_1': parameters['tau_1'],
        'tau_2': parameters['tau_2'],
        'tau_3': parameters['tau_3'],
        'tau_3_long': parameters['tau_3_long'],
        'tau_4': parameters['tau_4'],
        }
    
    # create random number generators
    seeds_generator=randomgen(seed+10000) # add to the seed so that rng3 and 4 differ from rng1 and 2 in Simplicity class
    # random number generators for population
    rng3 = randomgen(seeds_generator.integers(0,10000)) # for intra-host model states update
    rng4 = randomgen(seeds_generator.integers(0,10000)) # for electing individuals|lineages when reactions happen
    rng5 = randomgen(seeds_generator.integers(0,10000)) # for mutation model
    rng6 = randomgen(seeds_generator.integers(0,10000)) # for synthetic sequencing data
    
    # create population
    pop = Population(pop_size, I_0, ih_model_parameters, rng3,rng4,rng5,rng6, NSR_long,
                     long_shedders_ratio, sequence_long_shedders)
    return pop
        

