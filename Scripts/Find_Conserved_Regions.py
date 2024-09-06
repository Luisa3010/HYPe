####### get the conserved regions of a HYP1 G.Pallida read from a file and write
# them to a new file


import numpy as np
from tqdm import tqdm
from grep_motifs_semiglobal import semiglobal_matrix

# Load sequences from file
path = '/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/heads4_top_allele1_reads.fasta'
ids = []
seqs = []
with open(path, 'r') as file:
    for line in file:
        if line.startswith('>'):
            ids.append(line.strip('\n'))
        else:
            seqs.append(line.strip('\n'))

c1s =[]
c2s = []

# Use each of the two top sequences from heads4
for seq in tqdm(seqs):

    # Get conserved regions isolated
    # DNA sequences that mark the beginning and end of the HVD
    start = "GAAAGTGGTAAAAGACCCGGGAGC"
    end = "CATAAACACGGAGGTTATGACGAG"

    # Get the semiglobal alignments to determine start and end of the HVD
    HVD_start_scores = semiglobal_matrix(start, seq)[-1,:]
    HVD_end_scores = semiglobal_matrix(end, seq)[-1,:]
    HVD_start_pos = np.argmax(HVD_start_scores)
    HVD_end_pos = np.argmax(HVD_end_scores) - len(end)

    # Extract the conserved regions
    c1 = seq[:HVD_start_pos]
    c2 = seq[HVD_end_pos:]

    c1s.append(c1)
    c2s.append(c2)

# Write C1 and C2 to files
path = '/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/heads4_top_allele1_conserved1.fasta'
with open(path, 'w') as file:
    for region, id in zip(c1s, ids):
        file.write(id + '\n')
        file.write(region + '\n')

path = '/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/heads4_top_allele1_conserved2.fasta'
with open(path, 'w') as file:
    for region, id in zip(c2s, ids):
        file.write(id + '\n')
        file.write(region + '\n')

