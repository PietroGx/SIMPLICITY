"""
Created on Tue Jun  6 13:13:14 2023

@author: pietro
"""
import simplicity.host                as h
import simplicity.evolution.mutations as evo
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
                 size,I_0,
                 rng3,rng4,rng5,rng6,
                 reservoir=100000):
        
        # counter for inf reactions
        self.count_infections = 0
        self.count_infections_from_long_shedders = 0
        
        # random number generator
        self.rng3 = rng3
        self.rng4 = rng4
        self.rng5 = rng5
        self.rng6 = rng6 # for synthetic sequencing data
        
        # compartments and population groups
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
        
        # set attributes for evolutionary model
        
        self.ref_genome      = ref.get_reference()  # sequence of reference genome
        self.L_var           = len(self.ref_genome) # lenght of reference genome
        self.M               = evo.sub_matrix()     # substitution matrix
        
        self.active_variants_n = I_0                # number of IH viruses in the 
                                                    # infected population
        
        self.mutationsnumber = 0                    # counter of how many 
                                                    # substitutions happen in  
                                                    # one simulation run
        
        # track data for building the phylogenetic tree
        self.track_phylogeny = []
        self.phylodots = []
        self.lineage_frequency = {} # count lineages in the population
        
        self.lineages = {'wt': []} # dic of lineages sequences
        
        self.consensus_snapshot = [[[],1,0]] # [sequence, lineage_frequency, w_t(t_sim)]
        
        # system trajectory
        self.time         = 0 # time
        self.trajectory = [[self.time,
                           self.infected,
                           self.diagnosed,
                           self.recovered,
                           self.deceased,
                           self.infectious_normal,
                           self.detectables
                           ]]
        
        self.fitness_trajectory = [[0,[0,0]]]  
        self.sequencing_data = []
        
        # individuals
        self.individuals = {}              # store individuals data
        self.host_model  = h.Host()        # intra host model for normal individuals 
        
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
                     'viral_genomes'       : [[]],
                     'state_t'     : 0,
                     'state'       : 'susceptible',
                     'model'       : self.host_model,
                     'type'        : 'normal',
                     'fitness'     : 1,
                     'IH_virus_fitness' : [1],
                     'IH_virus_number'    : 1,
                     'IH_virus_names'   : [],
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
            dic[i]['IH_virus_names'] = ['wt']
            
            self.infected_i.append(self.susceptibles_i.pop(self.susceptibles_i.index(i)))
        # update variant-individuals track
        self.lineage_frequency[0] = {'wt':I_0}
        
        # return dictionary containing all individuals data (self.individuals)
        return dic
    
    def count_variants_t(self,t):
        '''
        Count how many individuals are infected by each variant at time
        t in the population

        Parameters
        ----------
        t : float
            time.
        '''
        # store variants here
        # {'variant name': number of infected patients}
        count_t = {}
        
        
        # once every 5 days save snapshot of lineages for consensus
        if np.floor(t%5) == 0:
            
            # Loop through the infected individuals and get a list of IH_virus names
            for key in self.infected_i:
                string_list = list(set(self.individuals[key]['IH_virus_names']))
            
                # Count each lineage and update count_t
                for s in string_list:
                    count_t[s] = count_t.get(s, 0) + 1
                    
            
            for key in count_t:
                
                self.consensus_snapshot.append([self.lineages[key],count_t[key
                                                   ]/sum(count_t.values()),t])
            
            self.lineage_frequency[t] = count_t
            
        # else just count lineages frequency    
        else: 
            # Loop through the infected individuals and get a list of IH_virus names
            for key in self.infected_i:
                string_list = list(set(self.individuals[key]['IH_virus_names']))
            
                # Count each lineage and update count_t
                for s in string_list:
                    count_t[s] = count_t.get(s, 0) + 1
            
            self.lineage_frequency[t] = count_t
        
    def diagnosis(self, seq_rate=0):
        '''
        Select an infected (and detectable) individual at random and tags it
        as "diagnosed". Update the infected and diagnosed compartments. 
        Update detectable_i, infectious_i, infected_i and diagnosed_i.
        '''
        self.infected  -= 1
        
        # select random patient to be diagnosed
        patient = self.rng4.choice(self.detectable_i)
        
        if self.rng6.uniform(0,1) < seq_rate:
            # store sequencing data
            # patient_id, time, genome, subst number, patient time, inf duration
            i = 0
            for genome in self.individuals[patient]['viral_genomes']:
                self.sequencing_data.append([
                    patient,
                    self.time,
                    genome,
                    len(genome),
                    self.individuals[patient]['type'],
                    self.time - self.individuals[patient]['t_infection'],
                    i
                    ])
                i+=1
        
        # update the active variants number
        self.active_variants_n -= self.individuals[patient]['IH_virus_number']
        
        # 0.01 of patients diagnosed die
        if self.rng4.uniform(0,1) < 0.01:
            self.deceased += 1
            
            # set patient as deceased
            self.individuals[patient]['t_final'] = self.time
            self.individuals[patient]['state']   = 'deceased'
            
            # remove individual's index from detectable_i, infectious_i and 
            # infected_i and add it to deceased_i
            try:
                ## in case the diagnosed patient was infectious, updated the infectious list
                self.infectious_normal_i.pop(self.infectious_normal_i.index(patient))
                self.infectious_normal -= 1
            except:
                pass
            ## update detectables
            self.detectable_i.pop(self.detectable_i.index(patient))
            self.detectables -= 1
            
            ## update deceased and infected
            self.deceased_i.append(self.infected_i.pop(self.infected_i.index(patient)))
           
        else:
            self.diagnosed += 1
            # set patient as diagnosed
            self.individuals[patient]['t_final'] = self.time
            self.individuals[patient]['state']   = 'diagnosed'
           
            # remove individual's index from detectable_i, infectious_i and 
            # infected_i and add it to diagnosed_i
            try:
                ## in case the diagnosed patient was infectious, updated the infectious list
                self.infectious_normal_i.pop(self.infectious_normal_i.index(patient))
                self.infectious_normal -= 1
            except:
                pass
            ## update detectables
            self.detectable_i.pop(self.detectable_i.index(patient))
            self.detectables -= 1
            ## update diagnosed and infected
            self.diagnosed_i.append(self.infected_i.pop(self.infected_i.index(patient)))
        
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
        # virus (select variant at random to be transmitted)
        index = self.rng4.integers(0,self.individuals[parent]['IH_virus_number'])
        self.individuals[new_inf]['viral_genomes'] = [copy.deepcopy(
            self.individuals[parent]['viral_genomes'][index])]
        # variant name
        self.individuals[new_inf]['IH_virus_names'].append(
            self.individuals[parent]['IH_virus_names'][index])
        # variant fitness
        self.individuals[new_inf]['IH_virus_fitness'] = [int(
            self.individuals[parent]['IH_virus_fitness'][index])]
        # update individual fitness (average of variants fitness)
        self.individuals[new_inf]['fitness'] = np.average(
                    self.individuals[new_inf]['IH_virus_fitness']) 
        # state
        self.individuals[new_inf]['state'] = 'infected'
        
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
            
    def update_time(self,time):
        # update the time 
        self.time = time
        
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
                           self.detectables
                           ])
    
    def add_variant(self):
        
        # select individual that will get a new variant
        new_var = self.rng4.choice(self.infected_i)
        
        # if variants present are less then the max
        if self.individuals[new_var]['IH_virus_number'] < self.individuals[new_var][
                                                               'IH_virus_max']:     
                ## duplicate variant at random
                # select randomly index of variant to duplicate
                index = self.rng4.integers(0,self.individuals[new_var][
                                                                   'IH_virus_number'])
                # add copy of variant to individual
                self.individuals[new_var]['viral_genomes'].append(copy.deepcopy(list(
                    self.individuals[new_var]['viral_genomes'][index])))  
                # add variant name
                self.individuals[new_var]['IH_virus_names'].append(
                    self.individuals[new_var]['IH_virus_names'][index]) 
                # add variant fitness 
                self.individuals[new_var]['IH_virus_fitness'].append(
                    self.individuals[new_var]['IH_virus_fitness'][index])
                # update individual fitness (average of variants fitness)
                self.individuals[new_var]['fitness'] = np.average(
                    self.individuals[new_var]['IH_virus_fitness'])
                # update the variants count
                self.individuals[new_var]['IH_virus_number'] +=1 
                # update the active variants number
                self.active_variants_n += 1
                
        # if variants present are equal to the max           
        elif self.individuals[new_var]['IH_virus_number'] == self.individuals[new_var
             ]['IH_virus_max'] and self.individuals[new_var]['IH_virus_max']>1:
                    ## delete a variant at random
                     # select randomly index of variant to delete
                    random_index = self.rng4.integers(0, len(
                        self.individuals[new_var]['viral_genomes']))
                     # delete variant
                    self.individuals[new_var]['viral_genomes'].pop(random_index)
                     # delete name
                    self.individuals[new_var]['IH_virus_names'].pop(random_index)
                     # delete fitness value
                    self.individuals[new_var]['IH_virus_fitness'].pop(random_index)
                    ## duplicate an existing variant at random
                     # select randomly index of variant to duplicate
                    index = self.rng4.integers(0,len(
                        self.individuals[new_var]['viral_genomes']))
                     # duplicate variant and append it to the variants list 
                    self.individuals[new_var]['viral_genomes'].append(copy.deepcopy(
                        list(self.individuals[new_var]['viral_genomes'][index])))    
                     # add variant name    
                    self.individuals[new_var]['IH_virus_names'].append(
                        self.individuals[new_var]['IH_virus_names'][index])
                     # add variant fitness 
                    self.individuals[new_var]['IH_virus_fitness'].append(
                        self.individuals[new_var]['IH_virus_fitness'][index])
                    # update individual fitness (average of variants fitness)
                    self.individuals[new_var]['fitness'] = np.average(
                        self.individuals[new_var]['IH_virus_fitness'])
        
        # raise error if there are too many variants
        elif self.individuals[new_var]['IH_virus_number'] > self.individuals[new_var]['IH_virus_max']:
                    raise ValueError("Too many variants in", new_var )
        else: 
            pass
            
    def mutate(self, e, dt, phenotype_model, consensus):
        '''
        Mutation model, mutates the viruses in the population.

        Parameters
        ----------
        e : float
            population evolutionary rate.
        dt : float
            delta t - extrande time step (in years)

        '''
        # select number of substitutions and positions in pooled genome
        select_pos = evo.select_positions(self, self.L_var, self.rng5, e, dt) 
        
        # if mutations are happening
        if select_pos != 'No substitutions':
            
            # positions in pooled genome
            positions = select_pos[0] 
            # vector of active variants
            active_variants = select_pos[1]
            # number of substitutions happening
            subst_number = select_pos[2]
            # count the number of substitutions happening during the simulation
            self.mutationsnumber += subst_number
            
            # map substitutions to index of variants
            subst_coord = evo.map_to_dic(active_variants, positions, self.L_var)
            
            # fetch the bases that will undergo substituion
            subst_coord = evo.fetch_bases(self, subst_coord)
            # count number of each nitrogenous base to be mutated (for bulk update)
            bases_count = evo.get_subst_numbers(subst_coord, subst_number)
            # use substitution matrix to select mutations that will happen
            unassigned_subst = evo.substitution(self.M, bases_count,self.rng5)
            # assign the mutations to their relative genome
            subst_coord = evo.assign_sub(unassigned_subst, subst_coord)
            # update variants in simulation with the corresponding substitution
            evo.update_variants(self, subst_coord)
            
            # update fitness score of individuals and variants
            # print('--X----X----X----X----X----X--')
            # print('Substitution coordinates: ', subst_coord)
            individuals_to_update = pheno.get_individuals_to_update(subst_coord)
            # print('Individuals to be updated: ', individuals_to_update)
            update_fitness = pheno.update_fitness_factory(phenotype_model)
            self.individuals = update_fitness(self.individuals, 
                                              individuals_to_update, consensus)
            # print('Mutation #', self.mutationsnumber)
            # print('')
        # else skip other steps
        
    def calculate_mean_fitness(self):
        value = []
        
        for i in self.infected_i:
                value.append(self.individuals[i]['fitness'])
                
        return [np.mean(value),np.std(value)]

    # store avg fitness
    def track_fitness(self):
        self.fitness_trajectory.append([self.time,self.calculate_mean_fitness()])
        
    def data(self):
        # return population dictionary as data frame
        df = pd.DataFrame(self.individuals).transpose()
        filtered_df = df[~df['state'].str.contains('susceptible')]
        return filtered_df
    
    def phylogeny(self):
        # return phylogeny dictionary as data frame
        return pd.DataFrame(self.track_phylogeny)
    
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
#####################################################################################################
    
def create_population(parameters):
    '''
    Create population instance from parameters file and return it.
    '''
    # population parameters
    pop_size = parameters['population_size']
    I_0      = parameters['infected_individuals_at_start']
    seed     = parameters['seed']
    
    # create random number generators
    seeds_generator=randomgen(seed+10000) # add to the seed so that rng3 and 4 differ from rng1 and 2 in Simplicity class
    # random number generators for population
    rng3 = randomgen(seeds_generator.integers(0,10000)) # for intra-host model states update
    rng4 = randomgen(seeds_generator.integers(0,10000)) # for electing individuals|variants when reactions happen
    rng5 = randomgen(seeds_generator.integers(0,10000)) # for mutation model
    rng6 = randomgen(seeds_generator.integers(0,10000)) # for synthetic sequencing data
    
    # create population
    pop = Population(pop_size, I_0, rng3,rng4,rng5,rng6)
    return pop
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
