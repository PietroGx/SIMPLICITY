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
import simplicity.phenotype.consensus  as c
from   simplicity.random_gen          import randomgen
import pandas as pd
import numpy as np
import scipy.stats
import scipy.special

# Internal switch, not a simulation parameter. The fitness trajectory is a
# diagnostic: nothing in the model reads it, and of the columns it writes only
# Time, Mean and Std are read downstream (plots_manager, the archived figures).
# Off by default; flip this to record it again.
track_fitness_traj = False

class Population:
    '''
    The class defines a population for the SIMPLICITY simulations. It contains
    the data about every individual as well as their intra-host model.
    '''
    def __init__(self,
                 size,I_0,
                 ih_model_parameters,
                 rng3,rng4,rng5,
                 NSR_long,
                 long_shedders_ratio=0,
                 sequence_long_shedders=False,
                 susceptibility_long=1.0,
                 write_fasta=False,
                 reservoir=100000):
        
        # random number generator
        self.rng3 = rng3 # for intra-host model states update
        self.rng4 = rng4 # for electing individuals|lineages when reactions happen
        self.rng5 = rng5 # for mutation model
        
        self.size = size
        # Share of the POPULATION that sheds long, fixed per individual at
        # creation (see _init_individuals) -- not a per-infection draw.
        self.long_shedders_ratio = long_shedders_ratio
        # Relative weight on who RECEIVES an infection; 1.0 = long shedders
        # are infected in proportion to their share of the susceptible pool.
        self.susceptibility_long = susceptibility_long

        # compartments and population groups ---------------------------------
        
        self.susceptibles = size - I_0 # number of susceptible individuals
        self.infected     = I_0        # number of infected individuals
        self.diagnosed    = 0          # number of diagnosed individuals
        self.diagnosed_standard = 0    #
        self.diagnosed_long     = 0    # 
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
        self.individuals = {}             #  store individuals data
        self.reservoir_i          = set() # set of indices of individuals in the reservoir
        
        self.susceptibles_i       = set() # set of susceptible individuals indices  
        
        self.infected_i           = set() # set of infected individuals indices 
        self.long_shedder_i       = set() # set of infected long shedder individuals
        
        self.diagnosed_i          = set() # set of diagnosed individuals indices
        self.diagnosed_standard_i = set() #
        self.diagnosed_long_i     = set() # 
        
        self.recovered_i          = set() # set of recovered individuals indices
        
        
        
        self.infectious_standard_i = set()
        self.infectious_long_i     = set()
        self.infectious_i          = set() # set of infectious individuals indices
        
        self.detectables_standard_i = set()
        self.detectables_long_i     = set()
        self.detectables_i          = set() # set of detectable individuals indices 
        
        self.exclude_i      = set()    # set to store the newly infected individual (excludes it from states update at time of infection)
        
        # set attributes for evolutionary model ===============================
        self.ref_genome      = ref.get_reference()  # sequence of reference genome
        self.L_lin           = len(self.ref_genome) # lenght of reference genome
        self.active_lineages_n = I_0                # number of IH viruses in the 
                                                    # infected population
        self.NSR_long = NSR_long # long shedders substitution rate
        
        # Phylogenetic tree data ----------------------------------------------
        self.sequence_long_shedders = sequence_long_shedders
        # FASTA output is opt-in: nothing in the repo reads either fasta file
        self.write_fasta = write_fasta
        self.phylogenetic_data = [{  'Time_emergence'  : 0,
                                     'Lineage_name'    : 'wt',
                                     'Lineage_parent'  : None,
                                     'Genome'          : {},
                                     'Host_type'       : 'standard',
                                     'Total_infections': 0
                                 }]
        self._phylo_name_map = { row['Lineage_name']: row for row in self.phylogenetic_data }
        self.phylodots = []         # needed to name lineages
        # ---------------------------------------------------------------------
        
        self.lineage_frequency = [] # count lineage frequency in the population
        
        # -------------------------------------------------------------------------   

        # Running weighted consensus. The per-entry history used to be kept as
        # a list and rescanned on every rebuild; it grew without bound (millions
        # of genome dicts by day 1095 at N=5000) and nothing ever read it back,
        # since only consensus_sequences_t is written out. See
        # consensus.ConsensusAccumulator.
        self.consensus = c.ConsensusAccumulator()
        self.consensus_sequences_t = [[{},0]] # list to store consensus everytime is calculated in a simulation [consensus,t]
        
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
        self.ih_model_parameters = ih_model_parameters
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
        # Long-shedder status is a host trait fixed here, not drawn at
        # infection: long_shedders_ratio is the share of the POPULATION that
        # sheds long, independent of infection status. Exact count over the
        # initial `size`; i.i.d. over the reservoir, which tops up the
        # susceptible pool on every diagnosis and would otherwise dilute it.
        n_initial_long = int(round(size * self.long_shedders_ratio))
        initial_long = set(self.rng3.choice(size, size=n_initial_long,
                                            replace=False).tolist()) \
                       if n_initial_long else set()
        # create an entry in the dictionary for each individual (0 to number of
        # total individuals in the population (reservoir))
        for i in range(self.reservoir):
            is_long = ((i in initial_long) if i < size
                       else self.rng3.uniform() < self.long_shedders_ratio)

            dic[i] = {
                     't_infection' : None,
                     't_not_infected': None,

                     't_infectious': None,
                     't_not_infectious': None,

                     't_diagnosis' : None,

                     'type'        : 'long_shedder' if is_long else 'standard',
                     'state_t'     : 0,
                     't_next_state': None,
                     'state'       : 'susceptible',
                     
                     'parent'      : None,
                     'inherited_lineage': None,
                     'new_infections'  : [],
                     
                     'IH_lineages'   : [],
                     'IH_unique_lineages_number': 0, #1
                     'IH_lineages_number'    : 0,    #1
                     'IH_lineages_max': (self.rng3.integers(5,16) if is_long
                                         else self.rng3.integers(1,5)),
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
            
        # set individuals infected at the beginning of the simulation.
        # The starting cohort already carries its own type from creation; that
        # sets its intra-host model (tau_3 vs tau_3_long) and, via
        # long_shedder_i, its mutation rate.
        for i in range(I_0): # update data of individuals infected at the beginning of the simulation

            ind_type = dic[i]['type']
            dic[i]['parent']       = 'root'
            dic[i]['t_infection']  = 0
            dic[i]['t_infectious'] = None
            dic[i]['state_t']      = 0
            # sample next jump time from exp dist.
            state_t = dic[i]['state_t']
            rate = - self.host_model[ind_type].A[state_t][state_t]
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
            if ind_type == 'long_shedder':
                self.long_shedder_i.add(i)
                self.long_shedders += 1

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
       return self._phylo_name_map[lineage_name]['Genome']
    # -------------------------------------------------------------------------
    #                               Updates
    # -------------------------------------------------------------------------
    
    def refresh_unique_lineages(self, i):
        '''Recompute a host's distinct lineage count and carry the delta into
        active_lineages_n, which sums distinct lineages over hosts. Call it
        wherever IH_lineages changes.'''
        ind = self.individuals[i]
        n = len(set(ind['IH_lineages']))
        self.active_lineages_n += n - ind['IH_unique_lineages_number']
        ind['IH_unique_lineages_number'] = n

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
    
                self.active_lineages_n -= ind['IH_unique_lineages_number']
    
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
        if not track_fitness_traj:
            return

        fitness_scores = [self.individuals[i]['fitness_score'] for i in self.infected_i]
        if not fitness_scores:
            self.fitness_trajectory.append([self.time, 0, 0, 0])
            return

        scores = np.asarray(fitness_scores, dtype=float)
        mean_fitness = scores.mean()
        std_fitness = scores.std()

        # Shannon entropy of the normalised fitness. This is what
        # scipy.stats.entropy does internally -- normalise, entr, sum -- without
        # its dispatch wrapper, and without the redundant second normalisation
        # the previous code paid for by pre-dividing.
        total = scores.sum()
        entropy = (float(np.sum(scipy.special.entr(scores / total)))
                   if total > 0 else 0.0)
    
        self.fitness_trajectory.append({
                                        'Time': self.time,
                                        'Mean': mean_fitness,
                                        'Std': std_fitness,
                                        'Entropy': entropy
                                    })
    
    def update_lineage_frequency_t(self, t):
        '''
        Frequency of each lineage among the infected, undiagnosed population at
        time t. Every host -- long shedders included -- contributes total weight
        1, split equally over the DISTINCT lineages it carries, so frequencies
        sum to 1 and a host carrying many lineages does not outvote one carrying
        few. Feeds both lineage_frequency.csv and the running consensus.
        '''
        share_lineages_t = {}   # lineage -> summed per-host share
        count_lineages_t = {}   # lineage -> hosts carrying it
        hosts_at_t = 0

        for individual_index in self.infected_i:
            unique_lineages = set(self.individuals[individual_index]['IH_lineages'])
            if not unique_lineages:
                continue
            hosts_at_t += 1
            host_share = 1.0 / len(unique_lineages)
            for lineage_name in unique_lineages:
                share_lineages_t[lineage_name] = share_lineages_t.get(lineage_name, 0.0) + host_share
                count_lineages_t[lineage_name] = count_lineages_t.get(lineage_name, 0) + 1

        if not hosts_at_t:
            return

        # Lineages carrying byte-identical genomes contribute identically to
        # the consensus -- add_lineage duplicates a lineage into a new slot and
        # the copy stays identical until it mutates -- and frequency enters
        # every weighted sum linearly, so they are merged into one accumulate
        # instead of one per lineage.
        by_genome = {}
        for lineage_name, share in share_lineages_t.items():
            frequency = share / hosts_at_t
            self.lineage_frequency.append({
                'Lineage_name'              : lineage_name,
                'Time_sampling'             : t,
                'Frequency_at_t'            : frequency,
                'Individuals_infected_at_t' : count_lineages_t[lineage_name],
                })
            genome = self.get_lineage_genome(lineage_name)
            key = tuple(sorted(genome.items()))
            if key in by_genome:
                by_genome[key][1] += frequency
            else:
                by_genome[key] = [genome, frequency]

        for genome, frequency in by_genome.values():
            self.consensus.add(genome, frequency, t)
            
# -----------------------------------------------------------------------------

    def update_ih_lineages_trajectories(self):
        # Only hosts that were ever infected carry a trajectory, and
        # individuals_data_to_df drops the rest, so normalise those alone. The
        # trajectory dicts live in self.individuals, so mutating them here is
        # exactly what the frame would have carried.
        kept = {i: ind for i, ind in self.individuals.items()
                if 'susceptible' not in ind['state']}
        for idx, row in kept.items():
            lineage_traj_dic = row['IH_lineages_trajectory']  
        
            for lineage in lineage_traj_dic:
                if lineage_traj_dic[lineage].get('ih_birth') is None:
                    lineage_traj_dic[lineage]['ih_birth'] = row['t_infection']
                if lineage_traj_dic[lineage].get('ih_death') is None:
                    lineage_traj_dic[lineage]['ih_death'] = row['t_not_infected']
                # normalize numpy scalars to plain floats for clean CSV
                # serialization; None left as-is (still-infected hosts
                # have no t_not_infected and keep ih_death = None)
                if lineage_traj_dic[lineage].get('ih_birth') is not None:
                    lineage_traj_dic[lineage]['ih_birth'] = float(lineage_traj_dic[lineage]['ih_birth'])
                if lineage_traj_dic[lineage].get('ih_death') is not None:
                    lineage_traj_dic[lineage]['ih_death'] = float(lineage_traj_dic[lineage]['ih_death'])
        
        return pd.DataFrame(kept).transpose()
        
    def individuals_data_to_df(self):
        # return population dictionary as data frame
        df = self.update_ih_lineages_trajectories()
        return df.drop('t_next_state', axis=1)
    
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
    sequence_long_shedders = parameters['sequence_long_shedders']
    susceptibility_long = parameters['susceptibility_long']
    write_fasta = parameters['write_fasta']
    
    NSR_long = parameters['nucleotide_substitution_rate_long']
    
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
    
    # create population
    pop = Population(pop_size, I_0, ih_model_parameters, rng3,rng4,rng5, NSR_long,
                     long_shedders_ratio, sequence_long_shedders,
                     susceptibility_long, write_fasta)
    return pop
        

