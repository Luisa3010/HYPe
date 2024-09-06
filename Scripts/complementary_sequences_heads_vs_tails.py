####### Take an equal (large) amount of sequences from heads and tails and look
# what sequences are present in heads but not in tails and vice versa to potentially find a
# difference between them


import os
import random
from Bio import SeqIO
from collections import Counter
import csv

# Function to sample sequences from a file
def sample_sequences(filepath, sample_size):
    with open(filepath, "r") as handle:
        sequences = list(SeqIO.parse(handle, "fasta"))
        if len(sequences) >= sample_size:
            sequences = random.sample(sequences, sample_size)

    sample = []
    for seq in sequences:
        cleaned_seq = seq.seq.replace('YQRGGG', 'YERGGG').replace('RDDQRG',
                                                               'RDNKRG').replace(
                'GNRGGG', 'SNRGGG').replace('SSRGD','SDRGD')
        sample.append(cleaned_seq)
        
    return sample




# Function to process directories
def process_directories(dir, head_output_file,
                        tail_output_file, sample_size = 3996):


    files = sorted(os.listdir(dir))
    head_files = [f'heads{i}.phred15.400to1400bp.alignscore2000.fulllength_aa_diff0' for i in range(1,11)]
    tail_files = [f'tails{i}.phred15.400to1400bp.alignscore2000.fulllength_aa_diff0' for i in range(1,11)]

    head_seqs = {}
    tail_seqs = {}

    for head_file, tail_file in zip(head_files, tail_files):
        head_path = os.path.join(dir, head_file)
        tail_path = os.path.join(dir, tail_file)

        
        heads_sample = sample_sequences(head_path, sample_size)
        tails_sample = sample_sequences(tail_path, sample_size)

        # Get all the sequences from heads and tails with their counts
        for seq in heads_sample:
            if not seq in head_seqs.keys():
                head_seqs[seq] = 1
            else:
                head_seqs[seq] += 1

        for seq in tails_sample:
            if not seq in tail_seqs.keys():
                tail_seqs[seq] = 1
            else:
                tail_seqs[seq] += 1

    # Calculate the complement counts
    for key in head_seqs.keys():
        if key in tail_seqs.keys():

            num_shared = min(head_seqs[key], tail_seqs[key])
            head_seqs[key] -= num_shared
            tail_seqs[key] -= num_shared


    # Write the complement to new files
    with open(head_output_file,'w' ) as outfile:
        outfile.write(f"Sequence,Count\n")
        for key, value in head_seqs.items():
            outfile.write(f"{key},{value}\n")


    with open(tail_output_file,'w' ) as outfile:
        outfile.write(f"Sequence,Count\n")
        for key, value in tail_seqs.items():
            outfile.write(f"{key},{value}\n")


# Define directories
output_heads = '/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/heads_complementary_seqs'
output_tails = '/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/tails_complementary_seqs'
dir = '/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_semiglobal'

# Run the process
process_directories(dir, output_heads, output_tails)









