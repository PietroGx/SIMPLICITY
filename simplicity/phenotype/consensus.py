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
    unique_positions = sorted({num for entry in data for num in entry[0]})
    num_index = {num: i for i, num in enumerate(unique_positions)}
    
    bases = ['A', 'T', 'C', 'G']  
    letter_index = {letter: i for i, letter in enumerate(bases)}
    
    # Initialize the matrix
    matrix = np.zeros((len(bases), len(unique_positions)), dtype=float)
    if not unique_positions:
        return matrix, bases, unique_positions

    # Every genome votes in every column: its own base where it stores one, the
    # reference otherwise. Since a genome differs from the reference at only a
    # handful of the columns, seed every column with the full weight on the
    # reference base and then correct where a genome differs. Identical result,
    # inner loop over mutations rather than over all columns.
    total = sum(n_infected * weight for _, n_infected, weight in data)
    ref_rows = np.fromiter((letter_index[reference[p]] for p in unique_positions),
                           dtype=np.intp, count=len(unique_positions))
    matrix[ref_rows, np.arange(len(unique_positions))] = total
    for seq, n_infected, weight in data:
        v = n_infected * weight
        for num, char in seq.items():
            col = num_index[num]
            matrix[ref_rows[col], col] -= v
            matrix[letter_index[char], col] += v
    
    return matrix, bases, unique_positions

def weighted_consensus(matrix, positions):
    '''
    calculate the weighted consensus sequence between individuals in the population
    '''
    bases = ['A', 'T', 'C', 'G']
    
    consensus = {}
    
    # Iterate through each column in the matrix to find the base with the highest count
    for col_idx, pos in enumerate(positions):
        # Use argmax to find the index of the maximum count in the column
        max_base_idx = np.argmax(matrix[:, col_idx])
        
        # Use the index to get the corresponding letter
        max_base= bases[max_base_idx]
        # if position is diffrent from wt, append it. (we encode sequences as only positions that differ from wt)
        if max_base != reference[pos]:
            # same encoding as a lineage genome: only positions differing from wt
            consensus[pos] = max_base
    
    return consensus

def consensus_distribution(matrix, bases, positions):
    '''Per-column distributions P(p,b), normalised column by column, and

        delta = sum over columns of [ P(p, consensus base) - P(p, reference) ]

    the offset used by the centred distributional distance.
    '''
    P = {}
    delta = 0.0
    for col, pos in enumerate(positions):
        total = matrix[:, col].sum()
        if not total:
            continue
        column = matrix[:, col] / total
        P[pos] = {b: float(column[i]) for i, b in enumerate(bases)}
        delta += float(column.max()) - P[pos][reference[pos]]
    return P, delta

def get_consensus(data,t):
    '''(consensus sequence, column distributions, delta). The argmax distance
    reads the first only; the distributional distance reads the other two.'''
    data = get_seq_weights(data,t)
    matrix, bases, positions= build_weighted_consensus_matrix(data)
    return (weighted_consensus(matrix, positions),
            *consensus_distribution(matrix, bases, positions))


# ## example use 
# data = [
#         [[[10,'A'],[45,'G']],10,10],
#         [[[13,'A']],       5,10],
#         [[[10,'C'],[45,'G']],10,10]
#         ]
# # data contains sequence, n_infected_t, t

# data = get_seq_weights(data,100)

# # Build the ndarray
# matrix, bases, positions= build_weighted_consensus_matrix(data)


# # Call the function with the matrix and unique numbers
# positions_with_max_letter = weighted_consensus(matrix, positions)


# # # To display the result similarly to a pandas DataFrame for visualization:
# # import pandas as pd
# # df = pd.DataFrame(matrix, index=bases, columns=positions)
# # print(df)

# # Display the result
# print(positions_with_max_letter)

# import simplicity.evolution.decoder as dec 

# example_consensus_decoded = dec.decode_genome(positions_with_max_letter)

# print(example_consensus_decoded)
