#!/usr/bin/env python

import os

# Define the input and output file paths
input_file = '/path/to/reference/mature.fa'
output_file = '/path/to/directory/mature_miRNA_sequence_dna.fa'

# Function to convert RNA sequences to DNA sequences, filtering only mmu- references
def convert_rna_to_dna(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                if 'mmu-' in line:
                    outfile.write(line)
            else:
                if 'mmu-' in previous_line:
                    dna_sequence = line.replace('U', 'T')
                    outfile.write(dna_sequence)
            previous_line = line

# Run the conversion
convert_rna_to_dna(input_file, output_file)

print("Conversion complete. The DNA sequences are saved in 'mature_miRNA_sequence_dna.fa'")


