##### get 20 reads that appear only once in the parsed data

from Bio import SeqIO
import random


# fasta_file = '/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT/deduced motifs/head1and10_and_tail1and10/heads1.phred15.400to1400bp.alignscore2000.fulllength_.AA'
# fasta_file = '/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT/fasta/heads1/aa_heads1_1'
# fasta_file = '/home/luisa/Documents/Work/Uni_Cambridge/Data/new_script_plasmid_mix1_rep1_dil0in1000_aa'
# fasta_file = '/home/luisa/Documents/Work/Uni_Cambridge/Data/new_script_plasmid_mix1_rep1_dil128in1000_aa'
# fasta_file = '/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_SW/HYP1_ONT_fasta/plasmid_mix1_rep1_dil0in1000.phred15.400to1400bp.alignscore2000.fulllength_dna'
# fasta_file = "/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/semiglobal_test_run/plasmid_mix1_rep1_dil0in1000.phred15.400to1400bp.alignscore2000.fulllength_aa"
fasta_file = "/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/semiglobal_test_run/plasmid_mix1_rep1_dil0in1000.phred15.400to1400bp.alignscore2000.fulllength_aa_n2843"


num_unique=20
sequence_count = {}
sequence_id = {}

# Parse the FASTA file and count occurrences of each sequence
for record in SeqIO.parse(fasta_file, "fasta"):
    sequence = str(record.seq)
    if sequence in sequence_count:
        sequence_count[sequence] += 1
    else:
        sequence_count[sequence] = 1
        sequence_id[sequence] = record.id

# Filter sequences that occur only once
unique_sequences = {}
for seq, count in sequence_count.items():
    if count == 1 or count == 2:
        unique_sequences[sequence_id[seq]] = seq

# Select 20 unique IDs randomly if there are more than 20

if len(unique_sequences) > num_unique:
    ids = random.sample(unique_sequences.keys(), num_unique)
else:
    ids = unique_sequences.keys()

# fasta_file2 = "/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT/fasta/heads1.phred15.400to1400bp.alignscore2000.fulllength.fasta"
# fasta_file2 = '/home/luisa/Documents/Work/Uni_Cambridge/Data/plasmids/HYP1_ONT_plasmid_fasta/plasmid_mix1_rep1_dil0in1000.phred15.400to1400bp.alignscore2000.fulllength.fasta'
# fasta_file2 = '/home/luisa/Documents/Work/Uni_Cambridge/Data/plasmids/HYP1_ONT_plasmid_fasta/plasmid_mix1_rep1_dil128in1000.phred15.400to1400bp.alignscore2000.fulllength.fasta'
# fasta_file2 = '/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_fasta/plasmid_mix1_rep1_dil0in1000.phred15.400to1400bp.alignscore2000.fulllength.fasta'
fasta_file2 = '/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/semiglobal_test_run/plasmid_mix1_rep1_dil0in1000.phred15.400to1400bp.alignscore2000.fulllength.fasta'

with open(fasta_file2, 'r') as file:
    is_in = False
    for line in file:
        if is_in:
            print(line.strip())
            print(unique_sequences[id])
            is_in = False
        if line.strip()[1:] in ids:
            is_in = True
            id = line.strip()[1:]
            print(id)


