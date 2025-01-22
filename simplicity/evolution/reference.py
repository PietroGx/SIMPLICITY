#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Pietro Gerletti
"""
import os
import simplicity.config as config

def get_reference():
    
    reference_filepath = os.path.join(config.get_data_dir(),'reference.txt')
    
    if not os.path.isfile(reference_filepath):
        from Bio import Entrez, SeqIO
        # Email address to identify yourself to the NCBI servers
        Entrez.email = "Gerlettip@rki.de"
        
        # Accession number of the sequence you want to fetch
        accession_number = "NC_045512.2"
        
        # Fetch the sequence using the accession number
        print("----------------------")
        print("----------------------")
        print('Importing SARS-CoV-2 reference genome from NCBI...')
        handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        print("----------------------")
        print("Genome imported:", record.description)
        spike_sequence = str(record.seq[20989:25956]) # fetch the spike gene
        handle.close()
    
        with open(reference_filepath, 'w') as file:
            # Write the string to the file
            file.write(spike_sequence)
        
    with open(reference_filepath, 'r') as file:
        # Read the content of the file
        spike_sequence = file.read()
        return spike_sequence