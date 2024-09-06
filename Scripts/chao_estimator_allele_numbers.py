####### calculate a chao estimator as a lower bound for the actuall number of alleles
# for heads and tails

import os.path
from statistics import median

from pydistinct.stats_estimators import *


# Get all the files
files = []
for i in range(1,11):
    file = f'/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_semiglobal/grouped/heads{i}.phred15.400to1400bp.alignscore2000.fulllength_aa_diff0'
    files.append(file)
    file = f'/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_semiglobal/grouped/tails{i}.phred15.400to1400bp.alignscore2000.fulllength_aa_diff0'
    files.append(file)


# Get all the Sequences from all the files
seqs = []
for infile in files:
    with open(infile, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                seqs.append(line.strip('\n'))



print('Observed Number of Alleles:')
print(len(list(set(seqs))))
print('Chao Estimator:')
print(chao_estimator(seqs))

