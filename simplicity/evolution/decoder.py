import simplicity.evolution.reference as ref

def decode_genome(single_encoded_genome):
    decoded_genome = ref.get_reference() # fetch reference genome
    if not isinstance(single_encoded_genome, list):
        raise ValueError(f"Expected a list, but got {type(single_encoded_genome)} instead.")
    
    if single_encoded_genome == []:
        return decoded_genome
    else:
        for substitution_coordinate in single_encoded_genome:            
            mutation_location = substitution_coordinate[0]
            mutated_Nbase     = substitution_coordinate[1]
            next_location_rejoin = substitution_coordinate[0]+1
            decoded_genome = decoded_genome[:mutation_location] + mutated_Nbase + decoded_genome[next_location_rejoin:]
        return decoded_genome