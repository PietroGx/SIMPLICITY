# -*- coding: utf-8 -*-
import simplicity.phenotype.weight as w 
import simplicity.evolution.reference as ref 
import numpy as np

reference = ref.get_reference()

def get_seq_weights(data,t_sim):
    """
    preprocess lineage data to calculate weighted consensus. 
    calculate weight for sequence at time t_sim
    
    Parameters
    ----------
    data : list 
        list of list containing lineages data.
        data contains lineage_name, sequence, n_infected_t, t
    t_sim : float
        time of the simulation at which we evalute the weights.

    Returns
    -------
    seq_for_consensus : list
        [sequence, #_infected_t, w_t(t_sim)].

    """
    # get parameters for w_t(t_sim)
    params = w.w_t_params()

    k_e = params[0]
    k_a = params[1]
    
    seq_for_consensus = []
    # loop over data, compute weight and add entry to list
    for lst in data:
        seq_for_consensus.append([
            lst[0],
            lst[1],
            w.weights(lst[2],t_sim, k_e, k_a, w.t_max)
            ])
   
    return seq_for_consensus

def build_weighted_consensus_matrix(data):
    '''
    

    Parameters
    ----------
    data : TYPE
        DESCRIPTION.

    Returns
    -------
    matrix : TYPE
        DESCRIPTION.
    bases : TYPE
        DESCRIPTION.
    unique_positions : TYPE
        DESCRIPTION.

    '''
    # Extract unique positions
    unique_positions= sorted(set(num for entry in data for num, _ in entry[0]))
    num_index = {num: i for i, num in enumerate(unique_positions)}
    
    bases = ['A', 'T', 'C', 'G']  
    letter_index = {letter: i for i, letter in enumerate(bases)}
    
    # Initialize the matrix
    matrix = np.zeros((len(bases), len(unique_positions)), dtype=float)
    
    # Populate the matrix, accounting for missing positions
    for seq, n_infected, weight in data:
        # Determine missing positions
        existing_positions = set(num for num, _ in seq)
        missing_positions = set(unique_positions) - existing_positions
        
        # Update matrix for existing positions
        for num, char in seq:
            row = letter_index[char]
            col = num_index[num]
            matrix[row, col] += n_infected * weight
        
        # Update matrix for missing positions
        for num in missing_positions:
            char = reference[num]  # fetch base from reference genome
            row = letter_index[char]
            col = num_index[num]
            matrix[row, col] += n_infected * weight 
    
    return matrix, bases, unique_positions

def weighted_consensus(matrix, positions):
    '''
    calculate the weighted consensus sequence between individuals in the population
    '''
    bases = ['A', 'T', 'C', 'G']
    
    consensus = []
    
    # Iterate through each column in the matrix to find the base with the highest count
    for col_idx, pos in enumerate(positions):
        # Use argmax to find the index of the maximum count in the column
        max_base_idx = np.argmax(matrix[:, col_idx])
        
        # Use the index to get the corresponding letter
        max_base= bases[max_base_idx]
        # if position is diffrent from wt, append it. (we encode sequences as only positions that differ from wt)
        if max_base != reference[pos]:
            # Append the position (num) and the letter with the highest count
            consensus.append([pos, max_base])
    
    return consensus

def get_consensus(data,t):
    data = get_seq_weights(data,t)
    matrix, bases, positions= build_weighted_consensus_matrix(data)
    return weighted_consensus(matrix, positions)


## example use 
data = [
        [[[10,'A'],[45,'G']],10,10],
        [[[13,'A']],       5,10],
        [[[10,'C'],[45,'G']],10,10]
        ]
# data contains sequence, n_infected_t, t

data = get_seq_weights(data,100)

# Build the ndarray
matrix, bases, positions= build_weighted_consensus_matrix(data)


# Call the function with the matrix and unique numbers
positions_with_max_letter = weighted_consensus(matrix, positions)


# # To display the result similarly to a pandas DataFrame for visualization:
# import pandas as pd
# df = pd.DataFrame(matrix, index=bases, columns=positions)
# print(df)

# Display the result
print(positions_with_max_letter)

import simplicity.evolution.decoder as dec 

example_consensus_decoded = dec.decode_genome(positions_with_max_letter)

print(example_consensus_decoded)
