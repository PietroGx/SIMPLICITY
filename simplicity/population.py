"""
Created on Tue Jun  6 13:13:14 2023

@author: pietro
"""
import simplicity.intra_host_model    as h
import simplicity.evolution.reference as ref
from   simplicity.random          import randomgen
import pandas as pd
import numpy as np

class Population:
    '''
    The class defines a population for the SIMPLICITY simulations. It contains
    the data about every individual as well as their intra-host model.
    '''
    def __init__(self,
                 size,I_0, tau_3,
                 rng3,rng4,rng5,rng6,
                 reservoir=100000):
        
        # counter for inf reactions
        self.count_infections = 0
        # self.count_infections_from_long_shedders = 0
        
        # random number generator
        self.rng3 = rng3 # for intra-host model states update
        self.rng4 = rng4 # for electing individuals|variants when reactions happen
        self.rng5 = rng5 # for mutation model
        self.rng6 = rng6 # for synthetic sequencing data
        
        # compartments and population groups ---------------------------------
        self.susceptibles = size-I_0 # number of susceptible individuals
        self.infected     = I_0 # compartment - number of infected individuals
        self.diagnosed    = 0   # compartment - number of diagnosed individuals
        self.recovered    = 0   # compartment - number of recovered individuals
        self.deceased     = 0   # compartment - number of deceased individuals
         
        self.reservoir = reservoir   # size of total population (not everyone is 
                                     # susceptible at the beginning, when 
                                     # individuals get removed from the system
                                     # new ones from the reservoir become 
                                     # susceptible) 
        
        self.infectious_normal   = 0        # number of infectious individuals
        self.detectables  = 0        # number of detectable individuals
        
        # set attributes for evolutionary model ===============================
        
        self.ref_genome      = ref.get_reference()  # sequence of reference genome
        self.L_lin           = len(self.ref_genome) # lenght of reference genome
        self.active_variants_n = I_0                # number of IH viruses in the 
                                                    # infected population
        
        # Phylogenetic tree data ----------------------------------------------
        self.phylogenetic_data = [{  'Time_emergence'  : 0,
                                     'Lineage_name'    : 'wt',
                                     'Lineage_parent'  : None,
                                     'Genome'          : [],
                                     'Host_type'       : 'normal'
                                     
                                 }]
        self.phylodots = []         # needed to name lineages
        self.sequencing_data = []
        # ---------------------------------------------------------------------
        
        self.lineage_frequency = [] # count lineage frequency in the population
        
        # =====================================================================
        
        self.consensus_snapshot = [[[],1,0]] # [sequence, lineage_frequency, w_t(t_sim)]
        self.consensus_sequences_t = [[[],0]] # list to store consensus everytime is calculated in a simulation [consensus,t]
        
        # system trajectory ---------------------------------------------------
        self.time         = 0 
        self.trajectory = [[self.time,
                           self.infected,
                           self.diagnosed,
                           self.recovered,
                           self.deceased,
                           self.infectious_normal,
                           self.detectables,
                           self.susceptibles
                           ]]
        
        self.last_infection = {}    # tracks the information about the last infection event
                                    # that happened in the simulaiton, used for 
                                    # R_effective calculations
        self.R_effective_trajectory = []
        self.fitness_trajectory = [[0,[0,0]]]  
        
        
        # individuals ---------------------------------------------------------
        self.individuals = {}              # store individuals data
        self.host_model  = h.Host(tau_3)   # intra host model for normal individuals 
        
        self.reservoir_i    = set()    # set of indices of individuals in the reservoir
        self.susceptibles_i = set()    # set of susceptible individuals indices  
        self.infected_i     = set()    # set of infected individuals indices 
        self.diagnosed_i    = set()    # set of diagnosed individuals indices
        self.recovered_i    = set()    # set of recovered individuals indices
        
        # self.infectious_i = []    
        self.infectious_normal_i = set()   # set of infectious individuals indices
        self.detectable_i        = set()   # set of detectable individuals indices 
        
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
                     't_infectious': None,
                     't_not_infectious': None,
                     
                     'parent'      : None,
                     'new_infections'  : [],
                     
                     'state_t'     : 0,
                     't_next_state': None,
                     'state'       : 'susceptible',
                     'model'       : self.host_model,
                     'type'        : 'normal',
                     
                     'IH_lineages'   : [],
                     'lineages_number': 1,
                     
                     'fitness'     : 1,
                     'IH_virus_fitness' : [1],
                     'IH_virus_number'    : 1,
                     'IH_virus_max': self.rng3.integers(1,6)
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
            rate = - dic[i]['model'].A[state_t][state_t]
            dic[i]['t_next_state'] = self.rng3.exponential(scale=1/rate) 
            
            dic[i]['state']        = 'infected'
            dic[i]['IH_lineages']  = ['wt']
            
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
        
    def get_next_ih_transition(self):
        ''' Get the next earliest time of an individual ih transition (for look ahead in extrande.)
        '''
        return min(
        (self.individuals[i]['t_next_state'] for i in self.infected_i if self.individuals[i]['t_next_state'] is not None),
        default=float('inf'))

    def update_states(self):
        """
        Update intra-host state for all infected individuals using the jump process.
        Handles:
          - advancing intra-host state
          - scheduling next transition
          - tracking infectious and detectable status
          - transitioning recovered individuals
        """
        
        to_update = [
            i for i in sorted(self.infected_i)
            if self.individuals[i]['t_next_state'] is not None and self.time >= self.individuals[i]['t_next_state']
        ]

        # Precompute transitions and rates
        ids = []
        scales = []
        
        for i in to_update:
            individual = self.individuals[i]
            individual['state_t'] += 1  # Move to next intra-host state
    
            # Update recovered individuals (state 20 = end of process) 
            if individual['state_t'] >= 20:
                individual['state'] = 'recovered'
                individual['t_not_infectious'] = self.time
                individual['t_next_state'] = None
    
                # Update compartments
                self.infected -= 1
                self.recovered += 1
                
                # Update sets
                self.infected_i.remove(i)
                self.recovered_i.add(i)
    
                # Remove from infectious/detectable sets if present
                self.infectious_normal_i.discard(i)
                self.detectable_i.discard(i)
    
                # Update active virus count
                self.active_variants_n -= individual['IH_virus_number']
    
                # Replace individual from reservoir
                new_susceptible_index = self.reservoir_i.pop()
                self.susceptibles_i.add(new_susceptible_index)
                self.susceptibles += 1
                continue  
    
            
            # Prepare for vectorized sampling of next intra-host transition time
            rate = -individual['model'].A[individual['state_t']][individual['state_t']]
            if rate > 0:
               ids.append(i)
               scales.append(1 / rate)
            else:
                raise ValueError(f"Invalid transition rate at state {individual['state_t']} for individual {i}")
    
            # Update infectious status (states 5–18) 
            if 4 < individual['state_t'] < 19 and individual['type'] == 'normal':
                if i not in self.infectious_normal_i:
                    self.infectious_normal_i.add(i)
                    individual['t_infectious'] = self.time
            else:
                if i in self.infectious_normal_i:
                    self.infectious_normal_i.remove(i)
                    individual['t_not_infectious'] = self.time
    
            # Update detectable status (states 5–19) 
            if 4 < individual['state_t'] < 20:
                self.detectable_i.add(i)
            else:
                self.detectable_i.discard(i)
        
        # Vectorized sampling of next intra-host transition times
        if ids:
            samples = self.rng3.exponential(scale=np.array(scales))
            for i, dt in zip(ids, samples):
                self.individuals[i]['t_next_state'] = self.time + dt
        
        # Update compartments
        self.infectious_normal = len(self.infectious_normal_i)
        self.detectables = len(self.detectable_i)

    def update_trajectory(self):
        # update the system trajectory
        self.trajectory.append([self.time,
                           self.infected,
                           self.diagnosed,
                           self.recovered,
                           self.infectious_normal,
                           self.detectables,
                           self.susceptibles
                           ])
    
    def update_fitness_trajectory(self):
        # update fitness trajectory
        def compute_mean_fitness(individuals, infected_i):
            value = []
            for i in infected_i:
                    value.append(individuals[i]['fitness'])
                    
            return [np.mean(value),np.std(value)]
     
        self.fitness_trajectory.append([self.time,compute_mean_fitness(self.individuals, self.infected_i)])
    
    def update_lineage_frequency_t(self, t):
        '''
        Count how many individuals are infected by each lineage at time
        t in the population and store the following information:
        time, lineage_name, frequency.
        '''
        # store lineages here: {'lineage name': number of infected individuals}
        count_lineages_t = {}
        
        # Loop through the infected individuals and get a list of unique IH_virus names (lineage)
        for individual_index in self.infected_i:
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
    # def update_R_effective_trajectory(self):
    #     for individual_index in self.infectious_normal_i:
    #         try:
    #             if self.time == self.last_infection['time_infection'] and individual_index == self.last_infection['transmitter_index']:
    #                 lineage_transmitted = self.last_infection['transmitted_lineage']
    #                 new_infections_at_t = 1
    #                 # store infection info for R effective
    #                 self.R_effective_trajectory.append([self.time,
    #                                            individual_index,
    #                                            lineage_transmitted,
    #                                            new_infections_at_t])
                
    #             else:
    #                 new_infections_at_t = 0
    #                 individual_lineages = self.individuals[individual_index]['IH_lineages']
    #                 # store infection info for R effective
    #                 for IH_lineage_index in range(0,len(individual_lineages)):
    #                     self.R_effective_trajectory.append([self.time,
    #                                                 individual_index,
    #                                                 individual_lineages[IH_lineage_index],
    #                                                 new_infections_at_t])
    #         except:
    #             if self.last_infection == {}:
    #                 new_infections_at_t = 0
    #                 individual_lineages = self.individuals[individual_index]['IH_lineages']
    #                 # store infection info for R effective
    #                 for IH_lineage_index in range(0,len(individual_lineages)):
    #                     self.R_effective_trajectory.append([self.time,
    #                                                 individual_index,
    #                                                 individual_lineages[IH_lineage_index],
    #                                                 new_infections_at_t])
                
    #             else:
    #                 raise ValueError('There is something wrong with R_effective_trajectory')
    # -------------------------------------------------------------------------
    
    def individuals_data_to_df(self):
        # return population dictionary as data frame
        df = pd.DataFrame(self.individuals).transpose()
        filtered_df = df[~df['state'].str.contains('susceptible')]
        return filtered_df
    
    def phylogenetic_data_to_df(self):
        # return phylogeny dictionary as data frame
        return pd.DataFrame(self.phylogenetic_data)
    
    def lineage_frequency_to_df(self):
        # return lineage_frequency as data frame
        return pd.DataFrame(self.lineage_frequency)
        
    def fitness_trajectory_to_df(self):
        times = [coord[0] for coord in self.fitness_trajectory]
        means = [coord[1][0] for coord in self.fitness_trajectory]
        stds = [coord[1][1] for coord in self.fitness_trajectory]
        
        # Create a pandas DataFrame
        df = pd.DataFrame({
            'Time': times,
            'Mean': means,
            'Std': stds
        })
        
        return df
    
    # def R_effective_trajectory_to_df(self):
    #     times = [item[0] for item in self.R_effective_trajectory]
    #     parents = [item[1] for item in self.R_effective_trajectory]
    #     lineages = [item[2] for item in self.R_effective_trajectory]
    #     new_infections_at_t = [item[3] for item in self.R_effective_trajectory]
        
    #     # Create a pandas DataFrame
    #     df = pd.DataFrame({
    #         'Time': times,
    #         'Individual': parents,
    #         'Lineage': lineages,
    #         'Infections_at_t':new_infections_at_t
    #     })
        
    #     return df

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
    tau_3    = parameters['tau_3']
    
    # create random number generators
    seeds_generator=randomgen(seed+10000) # add to the seed so that rng3 and 4 differ from rng1 and 2 in Simplicity class
    # random number generators for population
    rng3 = randomgen(seeds_generator.integers(0,10000)) # for intra-host model states update
    rng4 = randomgen(seeds_generator.integers(0,10000)) # for electing individuals|variants when reactions happen
    rng5 = randomgen(seeds_generator.integers(0,10000)) # for mutation model
    rng6 = randomgen(seeds_generator.integers(0,10000)) # for synthetic sequencing data
    
    # create population
    pop = Population(pop_size, I_0, tau_3, rng3,rng4,rng5,rng6)
    return pop
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
