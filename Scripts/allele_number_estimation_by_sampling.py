###### subsample the reads for heads and tails to get a curve to estimate the total
# number of possible alleles

import random
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
from scipy.optimize import curve_fit
import bisect
from tqdm import tqdm

def convert_global_to_local(file_paths, global_line_numbers):
    local_line_numbers = [[] for _ in file_paths]
    current_global_index = 0

    for i, file_path in enumerate(file_paths):
        with open(file_path, 'r') as file:
            num_lines = sum(1 for _ in file)

        local_lines_for_file = []
        for global_line_number in global_line_numbers:
            if current_global_index < global_line_number <= current_global_index + num_lines:
                local_lines_for_file.append(global_line_number - current_global_index)

        local_line_numbers[i] = local_lines_for_file
        current_global_index += num_lines

    return local_line_numbers


# read all the files ending in DNA from a folder
# folder_path = "/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT/deduced motifs/head1and10_and_tail1and10"
# dna_files = []
# for filename in os.listdir(folder_path):
#         if filename.endswith('DNA'):
#             file_path = os.path.join(folder_path, filename)
#             dna_files.append(file_path)

# dna_files = []
# for i in range(1,11):
#     dna_files.append(f"/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_SW/HYP1_ONT_fasta/heads{i}.phred15.400to1400bp.alignscore2000.fulllength_dna")
#
# dna_files = [f"/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_SW/HYP1_ONT_fasta/plasmid_mix1_rep1_dil0in1000.phred15.400to1400bp.alignscore2000.fulllength_dna",
# f"/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_SW/HYP1_ONT_fasta/plasmid_mix1_rep2_dil0in1000.phred15.400to1400bp.alignscore2000.fulllength_dna",
# f"/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_SW/HYP1_ONT_fasta/plasmid_mix1_rep3_dil0in1000.phred15.400to1400bp.alignscore2000.fulllength_dna"
#              ]


heads_files = []
for i in range(1,11):
    heads_files.append(f"/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_semiglobal/grouped/heads{i}.phred15.400to1400bp.alignscore2000.fulllength_aa_diff0")

# dna_files = []
# get number of sequences per file and the total number of sequences
num_lines = 0
file_sizes = []
for infile in heads_files:
    file_sizes.append(sum(1 for _ in open(infile)))

num_lines = sum(file_sizes)
num_seqs = round(num_lines/2)

# create different sample sizes
heads_sample_sizes = [i for i in range(1, num_seqs,10000)]
heads_sample_sizes.append(num_seqs -1)
# sample_sizes = [1,5,20,100,200,500,1000]

heads_num_uniques_per_sample = []
for sample_size in tqdm(heads_sample_sizes):
    # sample random sequences for the current sample size
    sample_ids = random.sample([i for i in range(1,num_seqs)], k = sample_size)
    sample_ids = sorted(sample_ids)

    for i in range(len(sample_ids)):
        sample_ids[i] = sample_ids[i] * 2

    # split sample ids by file
    local_sample_ids = convert_global_to_local(heads_files, sample_ids)

    last_file_end_id = 0
    unique_alleles = []
    for i, infile in enumerate(heads_files):
        with open(infile, 'r') as file:
            for j, line in enumerate(file) :
                if not line.startswith('>') and j + 1 in local_sample_ids[i]:
                    allele = line.strip().replace('YQRGGG','YERGGG').replace('RDDQRG', 'RDNKRG').replace('GNRGGG', 'SNRGGG').replace("SSRGD","SDRGD")
                    if not allele in unique_alleles:
                        unique_alleles.append(allele)

    heads_num_uniques_per_sample.append(len(unique_alleles))



tails_files = []
for i in range(1, 11):
    tails_files.append(
        f"/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_semiglobal/grouped/tails{i}.phred15.400to1400bp.alignscore2000.fulllength_aa_diff0")

# dna_files = []
# get number of sequences per file and the total number of sequences
num_lines = 0
file_sizes = []
for infile in tails_files:
    file_sizes.append(sum(1 for _ in open(infile)))

num_lines = sum(file_sizes)
num_seqs = round(num_lines / 2)

# create different sample sizes
tails_sample_sizes = [i for i in range(1, num_seqs, 10000)]
tails_sample_sizes.append(num_seqs - 1)
# sample_sizes = [1,5,20,100,200,500,1000]

tails_num_uniques_per_sample = []
for sample_size in tqdm(tails_sample_sizes):
    # sample random sequences for the current sample size
    sample_ids = random.sample([i for i in range(1, num_seqs)],
                               k = sample_size)
    sample_ids = sorted(sample_ids)

    for i in range(len(sample_ids)):
        sample_ids[i] = sample_ids[i] * 2

    # split sample ids by file
    local_sample_ids = convert_global_to_local(tails_files, sample_ids)

    last_file_end_id = 0
    unique_alleles = []
    for i, infile in enumerate(tails_files):
        with open(infile, 'r') as file:
            for j, line in enumerate(file):
                if not line.startswith('>') and j + 1 in local_sample_ids[
                    i]:
                    allele = line.strip().replace('YQRGGG',
                                                  'YERGGG').replace(
                        'RDDQRG', 'RDNKRG').replace('GNRGGG',
                                                    'SNRGGG').replace(
                        "SSRGD", "SDRGD")
                    if not allele in unique_alleles:
                        unique_alleles.append(allele)

    tails_num_uniques_per_sample.append(len(unique_alleles))


        # if the id of the line is in the list and the variant is not in the list yet, append it - hash list for quicker lookup?
        # with open(infile, 'r') as file:
        #     for i, line in enumerate(file):
        #         if not line.startswith('>') and i*2+1 in sample_ids:
        #             allele = line.strip()
        #             if not allele in unique_alleles:
        #                 unique_alleles.append(allele)
        #
        # num_uniques_per_sample.append(len(unique_alleles))
#
print(heads_sample_sizes)
print(heads_num_uniques_per_sample)
print(tails_sample_sizes)
print(tails_num_uniques_per_sample)
plt.plot(heads_sample_sizes, heads_num_uniques_per_sample)
plt.plot(tails_sample_sizes, tails_num_uniques_per_sample)
plt.legend(['heads','tails'])
plt.title('Allele variants for increasing sample sizes')
plt.xlabel('sample size')
plt.ylabel('number of unique variants')
plt.show()