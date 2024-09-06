##### split HVDs into their kmers to find the most common kmers - can be used to
# find potential motifs


def count_kmers(sequence, k, kmers):
    # Slide over the sequence and extract k-mers
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        if kmer in kmers:
            kmers[kmer] += 1
        else:
            kmers[kmer] = 1
    return kmers

def read_fasta(file_path):
    with open(file_path, 'r') as file:
        sequences = {}
        current_seq_id = None
        sequence = ''

        for line in file:
            line = line.strip()
            if line.startswith('>'):  # New sequence header
                if current_seq_id:
                    sequences[current_seq_id] = sequence
                # Extract the protein ID or some identifier from the header
                current_seq_id = line  # Takes the line as ID
                sequence = ''  # Reset sequence
            else:
                sequence += line.strip()  # Append sequence without newlines

        # Add the last sequence to the dictionary
        if current_seq_id:
            sequences[current_seq_id] = sequence

    return sequences


# read the data for the HYP protein sequences
#file_path = '/home/luisa/Documents/Work/Uni Cambridge/Data/HYP3_sanger/HYP3_protein_repeat_only_fasta (copy).txt'
#file_path = '/home/luisa/Documents/Work/Uni Cambridge/Data/HYP1_protein_repeat_only_fasta (copy).txt'
#file_path = '/home/luisa/Documents/Work/Uni Cambridge/Data/HYP1_sanger/HYP1_protein_repeat_only_fasta.txt'
file_path = '/home/luisa/Documents/Work/Uni Cambridge/Data/HYP3_sanger/HYP3_protein_repeat_only_fasta.txt'
family = "HYP-3"

sequences = read_fasta(file_path)


rf_lengths = range(4,30) # different kmer lengths
kmer_counts = {}

for seq_id, sequence in sequences.items():
    for k in rf_lengths:
        if family in seq_id:
            kmer_counts = count_kmers(sequence, k, kmer_counts)


# sort by most frequent letter combination
sorted_kmers = sorted(kmer_counts.items(), key=lambda item: item[1])

# print the results
for kmer, count in sorted_kmers:
    print(f"{kmer} {count}")


