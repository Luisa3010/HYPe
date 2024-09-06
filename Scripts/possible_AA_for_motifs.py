##### Get all the possible nucleotide combinations that could lead to the given
# amino acid sequenes to see what the motifs could possible be created from

import itertools
import pandas as pd

# Define the codon table
codon_table = {
    'A': ['GCT', 'GCC', 'GCA', 'GCG'],
    'C': ['TGT', 'TGC'],
    'D': ['GAT', 'GAC'],
    'E': ['GAA', 'GAG'],
    'F': ['TTT', 'TTC'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG'],
    'H': ['CAT', 'CAC'],
    'I': ['ATT', 'ATC', 'ATA'],
    'K': ['AAA', 'AAG'],
    'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'M': ['ATG'],
    'N': ['AAT', 'AAC'],
    'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'Q': ['CAA', 'CAG'],
    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    'W': ['TGG'],
    'Y': ['TAT', 'TAC'],
    '*': ['TAA', 'TAG', 'TGA']  # Stop codons
}

def generate_nucleotide_combinations(aa_sequence):
    # Get the list of codons for each amino acid in the sequence
    codon_lists = [codon_table[aa] for aa in aa_sequence]
    # Generate all possible combinations of codons
    all_combinations = [''.join(codon) for codon in itertools.product(*codon_lists)]
    return all_combinations

def find_variable_positions(combinations):
    variable_positions = [0]  * len(combinations[0])
    letter_variations = [[] for i in range(len(combinations[0]))]
    for combination in combinations:
        for i in range(0, len(combination)):
            if not combination[i] in letter_variations[i]:
                letter_variations[i].append(combination[i])
                variable_positions[i] += 1

    return variable_positions

def process_sequences(input_file, output_file):
    # Read the input sequences
    with open(input_file, 'r') as infile:
        sequences = [line.strip() for line in infile if line.strip()]

    # Prepare a list to store the results
    results = []

    # Generate nucleotide combinations and variable positions for each sequence
    for seq in sequences:
        combinations = generate_nucleotide_combinations(seq)
        variable_positions = find_variable_positions(combinations)
        for combination in combinations:
            results.append({'Amino Acid Sequence': seq,
                            'Nucleotide Combination': combination,
                            })
        results.append({'Amino Acid Sequence': seq,
                        'Nucleotide Combination': variable_positions
                        })

    # Convert the results to a DataFrame
    df = pd.DataFrame(results)

    # Save the DataFrame to a CSV file
    df.to_csv(output_file, index=False)


dir ='/home/luisa/Documents/Work/Uni Cambridge/Data/working_files/'
input_file = 'HYP3_motif_AA.txt'
output_file = 'HYP3_combinations.csv'

process_sequences(dir + input_file, dir + output_file)

print(f"Nucleotide combinations saved to {output_file}")
