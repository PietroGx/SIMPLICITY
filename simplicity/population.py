"""
Created on Tue Jun  6 13:13:14 2023

@author: pietro
"""
import simplicity.host                as h
import simplicity.evolution.mutations as evo_model
import simplicity.evolution.reference as ref
import simplicity.phenotype.update    as pheno
from   simplicity.random          import randomgen
import pandas as pd
import copy
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
        self.rng3 = rng3
        self.rng4 = rng4
        self.rng5 = rng5
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
        self.sub_matrix      = evo_model.sub_matrix()     # substitution matrix
        
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
        
        self.reservoir_i = []    # list of indices of individuals in the reservoir
        self.susceptibles_i = [] # list of susceptible individuals indices  
        self.infected_i = []     # list of infected individuals indices 
        self.diagnosed_i = []    # list of diagnosed individuals indices
        self.recovered_i = []    # list of recovered individuals indices
        self.deceased_i = []     # list of deceased individuals indices
        
        self.infectious_i = []   # list of infectious individuals indices 
        self.infectious_normal_i = []
        self.detectable_i = []   # list of detectable individuals indices 
        
        # dictionary with all individuals data 
        self.individuals = self._init_individuals(size,I_0)
        
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
                     't_final'     : None,
                     'parent'      : None,
                     'state_t'     : 0,
                     'state'       : 'susceptible',
                     'model'       : self.host_model,
                     'type'        : 'normal',
                     'IH_lineages'   : [],
                     
                     'fitness'     : 1,
                     'IH_virus_fitness' : [1],
                     'IH_virus_number'    : 1,
                     'IH_virus_max': self.rng3.integers(1,6),
                     'lineages_number': 1
                    }
            # add index of individuals to either susceptibles indices or to 
            # reservoir indices (the simulaiton starts with a pool of susceptibles)
            # that is replenished from the reservoir every time individuals
            # are removed from the system
            if i < size:
                self.susceptibles_i.append(i) 
            else:
                self.reservoir_i.append(i)
            
        # set individuals infected at the beginning of the simulation
        for i in range(I_0): # update data of individuals infected at the beginning of the simulation
            
            dic[i]['parent'] = 'root'
            dic[i]['t_infection'] = 0
            dic[i]['state_t'] = 5
            dic[i]['state'] = 'infected'
            dic[i]['IH_lineages'] = ['wt']
            
            self.infected_i.append(self.susceptibles_i.pop(self.susceptibles_i.index(i)))
        # update lineage_frequency
        self.lineage_frequency.append({'Lineage_name'              :'wt',
                                       'Time_sampling'             :0,
                                       'Frequency_at_t'            :1,
                                       'Individuals_infected_at_t' : I_0
                                       })
        
        # return dictionary containing all individuals data (self.individuals)
        return dic
    
    # -------------------------- SIR model reactions --------------------------
    
    def diagnosis(self, seq_rate=0):
        '''
        Select an infected (and detectable) individual at random and tags it
        as "diagnosed". Update the infected and diagnosed compartments. 
        Update detectable_i, infectious_i, infected_i and diagnosed_i.
        '''
        self.infected  -= 1
        
        # select random patient to be diagnosed
        patient_index = self.rng4.choice(self.detectable_i)
        
        if self.rng6.uniform(0,1) < seq_rate:
            # store sequencing data
            # patient_id, time, genome, subst number, patient time, inf duration
            i = 0
            for lineage_name in self.individuals[patient_index]['IH_lineages']:
                genome = self.get_lineage_genome(lineage_name)
                self.sequencing_data.append([
                    patient_index,
                    self.time,
                    genome,
                    len(genome),
                    self.individuals[patient_index]['type'],
                    self.time - self.individuals[patient_index]['t_infection'],
                    i
                    ])
                i+=1
        
        # update the active variants number
        self.active_variants_n -= self.individuals[patient_index]['IH_virus_number']
        
        # 0 % of patients diagnosed die
        if self.rng4.uniform(0,1) < 0:
            self.deceased += 1
            
            # set patient as deceased
            self.individuals[patient_index]['t_final'] = self.time
            self.individuals[patient_index]['state']   = 'deceased'
            
            # remove individual's index from detectable_i, infectious_i and 
            # infected_i and add it to deceased_i
            try:
                ## in case the diagnosed patient was infectious, updated the infectious list
                self.infectious_normal_i.pop(self.infectious_normal_i.index(patient_index))
                self.infectious_normal -= 1
            except:
                pass
            ## update detectables
            self.detectable_i.pop(self.detectable_i.index(patient_index))
            self.detectables -= 1
            
            ## update deceased and infected
            self.deceased_i.append(self.infected_i.pop(self.infected_i.index(patient_index)))
           
        else:
            self.diagnosed += 1
            # set patient as diagnosed
            self.individuals[patient_index]['t_final'] = self.time
            self.individuals[patient_index]['state']   = 'diagnosed'
           
            # remove individual's index from detectable_i, infectious_i and 
            # infected_i and add it to diagnosed_i
            try:
                ## in case the diagnosed patient was infectious, updated the infectious list
                self.infectious_normal_i.pop(self.infectious_normal_i.index(patient_index))
                self.infectious_normal -= 1
            except:
                pass
            ## update detectables
            self.detectable_i.pop(self.detectable_i.index(patient_index))
            self.detectables -= 1
            ## update diagnosed and infected
            self.diagnosed_i.append(self.infected_i.pop(self.infected_i.index(patient_index)))
        
        # add a susceptible back in from reservoir
        self.susceptibles += 1
        new_susceptible_index = self.reservoir_i.pop(0)
        self.susceptibles_i.append(new_susceptible_index)
       
    def infection(self):
        '''
        Select a random susceptible individual to be infected and tags it as 
        such. Update compartments and infected_i
        '''
        
        # select random patient to be the transmitter
        fitness_inf = [self.individuals[k]['fitness'] for k in self.infectious_normal_i]
        fitsum = np.sum(fitness_inf)
    
        if fitsum != 0:
            normed_fitness = [j/fitsum for j in fitness_inf]
            parent = self.rng4.choice(self.infectious_normal_i,p=normed_fitness)
        else:
            parent = self.rng4.choice(self.infectious_normal_i)
            
        # count inf reactions
        self.count_infections += 1
        
        # select random patient to be infected
        new_inf = self.rng4.choice(self.susceptibles_i)
            
        # update the active variants number
        self.active_variants_n += self.individuals[new_inf]['IH_virus_number']
        
        # move index to infected list
        self.infected_i.append(self.susceptibles_i.pop(
            self.susceptibles_i.index(new_inf)))
        # time of infection
        self.individuals[new_inf]['t_infection'] = self.time
        # parent
        self.individuals[new_inf]['parent'] = parent
        
        # virus (select lineage at random to be transmitted)
        index = self.rng4.integers(0,self.individuals[parent]['IH_virus_number'])
        transmitted_lineage = self.individuals[parent]['IH_lineages'][index]
        self.individuals[new_inf]['IH_lineages'].append(transmitted_lineage)
        # lineage fitness
        self.individuals[new_inf]['IH_virus_fitness'] = [
            self.individuals[parent]['IH_virus_fitness'][index]]
        
        # update individual fitness (average of variants fitness)
        self.individuals[new_inf]['fitness'] = np.average(
                    self.individuals[new_inf]['IH_virus_fitness']) 
        
        # state
        self.individuals[new_inf]['state'] = 'infected'
        
        # store infection info for R effective
        self.last_infection = {'time_infection':self.time, 
                               'transmitter_index':parent,
                               'transmitted_lineage':transmitted_lineage}
        
        # update susceptibles and infected 
        self.susceptibles -= 1
        self.infected     += 1
        
    def recovery(self):
        '''
        Tag recovered individuals as such, update the infected and recovered
        compartments and add the index of recovered individuals to recovered_i.
        '''
        for i in self.infected_i:
            if self.individuals[i]['state_t'] > 19:
                self.individuals[i]['t_final'] = self.time
                self.individuals[i]['state']   = 'recovered'
                
                # add index to recovered and remove it from infected indices
                self.recovered_i.append(self.infected_i.pop(self.infected_i.index(i)))
                # update the active variants number
                self.active_variants_n -= self.individuals[i]['IH_virus_number']
                    
                self.infected  -= 1 
                self.recovered += 1
                # add a susceptible back in from reservoir
                self.susceptibles += 1
                # add a new individual
                new_susceptible_index = self.reservoir_i.pop(0)
                self.susceptibles_i.append(new_susceptible_index)
                
    def add_variant(self):
        
        # select individual that will get a new variant
        individual_index = self.rng4.choice(self.infected_i)
        
        # if variants present are less then the max
        if self.individuals[individual_index]['IH_virus_number'] < self.individuals[individual_index][
                                                               'IH_virus_max']:     
                ## duplicate variant at random
                # select randomly index of lineage to duplicate
                index = self.rng4.integers(0,self.individuals[individual_index][
                                                                   'IH_virus_number'])
                # add IH lineage name
                self.individuals[individual_index]['IH_lineages'].append(
                    self.individuals[individual_index]['IH_lineages'][index]) 
                # add IH lineage fitness 
                self.individuals[individual_index]['IH_virus_fitness'].append(
                    self.individuals[individual_index]['IH_virus_fitness'][index])
                # update individual fitness (average of variants fitness)
                self.individuals[individual_index]['fitness'] = np.average(
                    self.individuals[individual_index]['IH_virus_fitness'])
                # update the variants count
                self.individuals[individual_index]['IH_virus_number'] +=1 
                # update the active variants number
                self.active_variants_n += 1
                
        # if variants present are equal to the max           
        elif self.individuals[individual_index]['IH_virus_number'] == self.individuals[individual_index
             ]['IH_virus_max'] and self.individuals[individual_index]['IH_virus_max']>1:
                    ## delete a variant at random
                     # select randomly index of variant to delete
                    random_index = self.rng4.integers(0, len(
                        self.individuals[individual_index]['IH_lineages']))
                     # delete lineage
                    self.individuals[individual_index]['IH_lineages'].pop(random_index)
                     # delete fitness value
                    self.individuals[individual_index]['IH_virus_fitness'].pop(random_index)
                    ## duplicate an existing IH lineage at random
                     # select randomly index of variant to duplicate
                    index = self.rng4.integers(0,len(
                        self.individuals[individual_index]['IH_lineages']))

                     # add lineage name    
                    self.individuals[individual_index]['IH_lineages'].append(
                        self.individuals[individual_index]['IH_lineages'][index])
                     # add IH lineage fitness 
                    self.individuals[individual_index]['IH_virus_fitness'].append(
                        self.individuals[individual_index]['IH_virus_fitness'][index])
                    # update individual fitness (average of IH lineages fitness)
                    self.individuals[individual_index]['fitness'] = np.average(
                        self.individuals[individual_index]['IH_virus_fitness'])
        
        # raise error if there are too many variants
        elif self.individuals[individual_index]['IH_virus_number'] > self.individuals[individual_index]['IH_virus_max']:
                    raise ValueError("Too many variants in", individual_index )
        else: 
            pass
 
    # -------------------------------------------------------------------------
    #                              Mutations
    # -------------------------------------------------------------------------
    
    def get_lineage_genome(self, lineage_name):
       '''
       Fetch lineage genome from lineage name
       '''
       genome = next((d['Genome'] for d in self.phylogenetic_data if d['Lineage_name'] == lineage_name))
       # print(f'{lineage_name} Genome: ', genome)
       return genome
   
    def mutate(self, e, dt, phenotype_model, *args):
        '''
        Mutation model, mutates the viruses in the population.

        Parameters
        ----------
        e : float
            population nucleotide substitution rate.
        dt : float
            delta t - extrande time step (in years)

        '''
        
        # select number of substitutions and positions in pooled genome
        select_pos = evo_model.select_positions(self, self.L_lin, self.rng5, e, dt) 
        
        # if mutations are happening
        if select_pos != 'No substitutions':
            
            # positions in pooled genome
            positions = select_pos[0] 
            # vector of active variants
            active_lineages = select_pos[1]
            # number of substitutions happening
            subst_number = select_pos[2]
            
            # map substitutions to index of variants
            subst_coord = evo_model.map_to_dic(active_lineages, positions, self.L_lin)
            
            # fetch the bases that will undergo substituion
            subst_coord = evo_model.fetch_bases(self, subst_coord)
            # count number of each nitrogenous base to be mutated (for bulk update)
            bases_count = evo_model.get_subst_numbers(subst_coord, subst_number)
            # use substitution matrix to select mutations that will happen
            unassigned_subst = evo_model.substitution(self.sub_matrix, bases_count,self.rng5)
            # assign the mutations to their relative genome
            subst_coord = evo_model.assign_sub(unassigned_subst, subst_coord)
            # update variants in simulation with the corresponding substitution
            evo_model.update_lineages(self, subst_coord)
            
            # update fitness score of individuals and variants
            
            # print('--X----X----X----X----X----X--')
            # print('Substitution coordinates: ', subst_coord)
            individuals_to_update = pheno.get_individuals_to_update(subst_coord)
            # print('Individuals to be updated: ', individuals_to_update)
            update_fitness = pheno.update_fitness_factory(phenotype_model)
            if args:
                consensus = args[0]
                update_fitness(self, individuals_to_update, consensus)
            else:
                update_fitness(self, individuals_to_update)
            # print('')
        # else skip other steps
    
    # -------------------------------------------------------------------------
    #                               Updates
    # -------------------------------------------------------------------------
    
    def update_time(self,time):
        # update the time 
        self.time = time
        
    def update_states(self,delta_t):
        '''
        Update intra-host model states for each individual in the simulation.
        '''
        # stores p_t for each state
        self.host_model.probabilities_t(delta_t)
        
        # draw random variables for each infected individual in the population
        tau = self.rng3.uniform(size=len(self.infected_i))
        
        i =0
        for key in self.infected_i:
            state = self.individuals[key]['state_t']
            self.individuals[key]['state_t'] = self.individuals[key]['model'
             ].update_state(self.individuals[key]['model'
             ].probabilities[state],state,tau[i])
            i += 1       
        
    def update_inf_det(self):
        '''
        Update the lists of infectious and detectables individuals.
        '''
        # remove indeces of individuals not infectious anymore (NORMAL)
        for i in self.infectious_normal_i:
            if self.individuals[i]['state_t']>=19:
                self.infectious_normal_i.pop(self.infectious_normal_i.index(i))
        # update the list of infectious individuals      
        for i in list(set(self.infected_i).symmetric_difference(set(self.infectious_normal_i))):
            if self.individuals[i]['state_t']>4 and self.individuals[i]['state_t']<19 and self.individuals[i]['type']=='normal':
                self.infectious_normal_i.append(i)
        
        # remove indices of individuals not detectable anymore
        for i in copy.deepcopy(self.detectable_i):
            if self.individuals[i]['state_t']==20:
                self.detectable_i.pop(self.detectable_i.index(i))
        # update the list of detectable individuals
        for i in list(set(self.infected_i).symmetric_difference(set(self.detectable_i))):
            if self.individuals[i]['state_t']>4 and self.individuals[i]['state_t']<20:
                self.detectable_i.append(i)

        self.detectables = len(self.detectable_i)
        self.infectious_normal  = len(self.infectious_normal_i)
        
    def update_trajectory(self):
        # update the system trajectory
        self.trajectory.append([self.time,
                           self.infected,
                           self.diagnosed,
                           self.recovered,
                           self.deceased,
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
    def update_R_effective_trajectory(self):
        for individual_index in self.infected_i:
            try:
                if self.time == self.last_infection['time_infection']:
                    individual_index = self.last_infection['transmitter_index']
                    lineage_transmitted = self.last_infection['transmitted_lineage']
                    new_infections_at_t = 1
                    # store infection info for R effective
                    self.R_effective_trajectory.append([self.time,
                                               individual_index,
                                               lineage_transmitted,
                                               new_infections_at_t])
                
                else:
                    new_infections_at_t = 0
                    individual_lineages = self.individuals[individual_index]['IH_lineages']
                    # store infection info for R effective
                    for IH_lineage_index in range(0,len(individual_lineages)):
                        self.R_effective_trajectory.append([self.time,
                                                    individual_index,
                                                    individual_lineages[IH_lineage_index],
                                                    new_infections_at_t])
            except:
                if self.last_infection == {}:
                    new_infections_at_t = 0
                    individual_lineages = self.individuals[individual_index]['IH_lineages']
                    # store infection info for R effective
                    for IH_lineage_index in range(0,len(individual_lineages)):
                        self.R_effective_trajectory.append([self.time,
                                                    individual_index,
                                                    individual_lineages[IH_lineage_index],
                                                    new_infections_at_t])
                
                else:
                    raise ValueError('There is something wrong with R_effective_trajectory')
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
    
    def R_effective_trajectory_to_df(self):
        times = [item[0] for item in self.R_effective_trajectory]
        parents = [item[1] for item in self.R_effective_trajectory]
        lineages = [item[2] for item in self.R_effective_trajectory]
        new_infections_at_t = [item[3] for item in self.R_effective_trajectory]
        
        # Create a pandas DataFrame
        df = pd.DataFrame({
            'Time': times,
            'Individual': parents,
            'Lineage': lineages,
            'Infections_at_t':new_infections_at_t
        })
        
        return df

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
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
