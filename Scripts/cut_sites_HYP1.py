###### create and store a pandas df with information about potential cut sites
# in rare alleles to use for further analysis


import pandas as pd
from collections import Counter
import os
from tqdm import tqdm

# Get all the file names
directory = '/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_semiglobal/grouped'
files = []
for i in range(1,11):
    path = f'tails{i}.phred15.400to1400bp.alignscore2000.fulllength_aa_diff0'
    files.append(os.path.join(directory,path))
    path = f'heads{i}.phred15.400to1400bp.alignscore2000.fulllength_aa_diff0'
    files.append(os.path.join(directory, path))


#
# files = []
# for file in os.listdir(directory):
#     if os.path.isfile(os.path.join(directory, file)) and 'aa' in file:
#         files.append(file)



column_names = [
        'Allele1',
        'Allele2',
        'TargetAllele',
        'Allele1CutPos',
        'Allele2CutPos',
        'Allele1MotifBefore',
        'Allele2MotifBefore',
        'Allele1MotifAfter',
        'Allele2MotifAfter',
        'Id',
        'Rep',
        'Category'
]
df = pd.DataFrame(columns = column_names)

for file in tqdm(files):

    # Read all the sequences with their ids from the files
    seqs = []
    ids = []
    with open(file, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                seqs.append(line.strip('\n'))
            else:
                ids.append(line.strip('>\n'))

    # get the counts of the alleles to determine the top most present alleles
    counts = Counter(seqs)
    top_alleles = []
    for seq, count in counts.most_common(10):
        if count / len(seqs) > 0.01:
            top_alleles.append(seq.split(' '))

    # check which combinations from cut up top alleles can be used to create each sequence
    for seq, id in zip(seqs,ids):
        target_seq = seq.split(' ')
        # only check if the sequence is not one of the top alleles
        if not target_seq in top_alleles:
            for allele1 in top_alleles:
                for allele2 in top_alleles:
                    # get the start and the ends of the two original sequences
                    for i in range(len(allele1)):
                        for j in range(len(allele2)):
                            start = allele1[:i]
                            end = allele2[j:]
                            combination = []
                            combination.extend(start)
                            combination.extend(end)
                            if combination == target_seq:
                                df = df._append({
                                        'Allele1': ' '.join(allele1), # top allele used for the start
                                        'Allele2': ' '.join(allele2), # top allele used for the end
                                        'TargetAllele': ' '.join(combination), # target sequence
                                        'Allele1CutPos': i,  # cut position allele1
                                        'Allele2CutPos': j,  # cut position allele2
                                        'Allele1MotifBefore': allele1[i - 1], # motif before the cut in allele1
                                        'Allele2MotifBefore': allele1[i], # motif before the cut in allele2
                                        'Allele1MotifAfter': allele2[j - 1], # motif after the cut in allele1
                                        'Allele2MotifAfter': allele2[j], # motif after the cut in allele2
                                        'Id': id, # Target sequence Id
                                        'Rep': os.path.basename(file.name), # file name
                                        'Category': 'heads' if 'heads' in os.path.basename(file.name) else 'tails' # heads ot tails
                                }, ignore_index = True)



print(df)

df.to_csv('/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/HYP1_HVD_cut_sites.csv')
