# -*- coding: utf-8 -*-
'''
Here there are the functions we use to calculate the hamming distance for the 
phenotype model. Please consider that they are adapted to the data sctructure we
use to store genomic data (we only store the positions and mutations that are 
                           different from the reference genome.)
'''
from simplicity.evolution import reference as ref
reference = ref.get_reference()

def hamming(lineage):
    # compute the hamming distance of a lineage from reference genome
    distance = 1
    
    for mutation in lineage:
        if reference[mutation[0]] != mutation[1]:
            distance += 1
    return distance

def hamming_true(lineage):
    # compute the hamming distance of a lineage from reference genome
    distance = 0
    
    for mutation in lineage:
        if reference[mutation[0]] != mutation[1]:
            distance += 1
    return distance

def hamming_iw(lineage,consensus):
    # compute the hamming distance of a lineage from consensus sequence
    
    # Convert sequences into dictionaries
    dict1 = {position: base for position, base in lineage}
    dict2 = {position: base for position, base in consensus}
    
    # Get the set of all unique positions in either sequence
    all_positions = set(dict1.keys()) | set(dict2.keys())
    
    # Initialize distance
    distance = 0
    
    # Iterate through all positions and compare bases
    for position in all_positions:
        base1 = dict1.get(position)
        base2 = dict2.get(position)
        
        # If bases are different or only one sequence contains the position, increment distance
        if base1 != base2:
            distance += 1
    
    return distance

    
