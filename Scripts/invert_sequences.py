###### Get the complementary inverted strands for a fasta file of DNA sequence


def reverse_complement(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reversed_sequence = sequence[::-1]  # Reverse the sequence
    reversed_complemented_sequence = ''.join(complement[base] for base in reversed_sequence)
    return reversed_complemented_sequence

def process_fasta(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        sequence = ''
        header = ''
        for line in infile:
            if line.startswith('>'):
                if sequence:
                    # Write the processed sequence to the output file
                    reversed_complement = reverse_complement(sequence)
                    outfile.write(f'{header}\n{reversed_complement}\n')
                header = line.strip()  # Save the header
                sequence = ''  # Reset the sequence
            else:
                sequence += line.strip()
        # Process and write the last sequence
        if sequence:
            reversed_complement = reverse_complement(sequence)
            outfile.write(f'{header}\n{reversed_complement}\n')

# Replace 'input.fasta' with your input file name and 'output.fasta' with your desired output file name
dir = "/home/luisa/Documents/Work/Uni Cambridge/Data/HYP1_ONT_sample/"
input_file = 'H110.align.fasta'
output_file = 'H110_align_reversed.fasta'

process_fasta(dir + input_file, dir + output_file)

print(f"Processed sequences saved to {output_file}")
