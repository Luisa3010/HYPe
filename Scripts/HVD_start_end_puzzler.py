###### for a given file of parsed HVDs, get the most frequent alleles and check
# if the remaning alleles can be puzzled together from starts and ends of those


from collections import Counter
import Levenshtein as lstn
import pandas as pd
import numpy as np
import os
from tqdm import tqdm


# path = '/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_semiglobal/grouped/plasmid_mix2_rep2_dil0in1000.phred15.400to1400bp.alignscore2000.fulllength_aa_diff0'
# path = '/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP_all_sanger/HYP1_sanger/HYP1_protein_repeat_only_cleaned_stripped.txt'

directory = '/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_semiglobal/grouped'
#
# files = []
# for file in os.listdir(directory):
#     if os.path.isfile(os.path.join(directory, file)) and 'aa' in file:
#         files.append(file)

files = []
for i in range(1,11):
    path = f'tails{i}.phred15.400to1400bp.alignscore2000.fulllength_aa_diff0'
    files.append(os.path.join(directory,path))
    path = f'heads{i}.phred15.400to1400bp.alignscore2000.fulllength_aa_diff0'
    files.append(os.path.join(directory, path))

col_labels = [ 'Reads', 'RareVariants', 'TemplateSwitches', 'SNPs', 'TemplateSwitch+SNPs', 'NotExplained' ]
results = pd.DataFrame(np.zeros((len(files),len(col_labels))), columns = col_labels, index = files)





for path in files:

    #
    n = 0
    with open(os.path.join(directory, path), 'r') as file:
        seqs = []
        for line in file:
            if not line.startswith('>'):
                seqs.append(line.strip('\n'))
                n +=1

    # Get top sequences in the files - over 1%
    top_seqs = []
    top_n = 0
    counts = Counter(seqs)
    for seq, count in counts.most_common(10):
        if count/n > 0.01:
            top_seqs.append(seq)
            top_n += count

    # known top_sequences in heads and tails
    # top_seqs.extend(['YERGGG SDRGGG RDRGD YERGGG SNRGGG SNRGGG SNRGGG SNRGGG RDRGD YERGGG SNRGGG SDRGD YERGGG SDRGGG RDRGD YERGGG SNRGGG SNRGGG REGGD YERGGG SNRGGG RDRGD YERGGG',
    #             'YERGGG SDRGGG RDNKRG SDRGD YERGGG SDRGGG SNRGGG SDRGD YERGGG SNRGGG SNRGGG SDRGD YERGGG SNRGGG SDRGD YERGGG SNRGGG SNRGGG SNRGGG REGGD YERGGG SNRGGG RDRGD YERGGG',
    #             'YERGGG RDNKRG SNRGGG RDRGD YERGGG SDRGGG SDRGGG RDRGD YERGGG SNRGGG SNRGGG RDRGD YERGGG SNRGGG SDRGD YERGGG SNRGGG SDRGD YERGGG SNRGGG RDRGGG REGGD YERGGG SNRGGG RDRGD YERGGG',
    #             'YERGGG SDRGGG RDRGD YERGGG SNRGGG SNRGGG RDRGD YERGGG SNRGGG SDRGD YERGGG SNRGGG SNRGGG SDRGD YERGGG SNRGGG SNRGGG REGGD YERGGG SNRGGG RDRGD YERGGG',
    #             'YERGGG SDRGGG RDNKRG SNRGGG SDRGD YERGGG SNRGGG REGGD YERGGG SNRGGG RDRGD YERGGG SNRGGG RDRGD YERGGG SNRGGG SDRGE YERGGG SNRGGG SDRGD YERGGG SNRGGG RDRGD YERGGG'
    #             ])


    # Get the starts and the ends of the top sequences as potential building blocks for HVDs
    starts = []
    ends = []
    for seq in top_seqs:
        # cut between any motifs

        motifs = seq.split(' ')
        for i in range(len(motifs)):
            starts.append((' ').join(motifs[:i]) + ' ')
            ends.append((' ').join(motifs[i:]))
    # cut after YERGGG
    #     blocks = seq.split('YERGGG ')[1:]
    #     for i in range(len(blocks)):
    #         starts.append('YERGGG ' + 'YERGGG '.join(blocks[:i]))
    #         ends.append('YERGGG ' + 'YERGGG '.join(blocks[i:]))
    ends.append('YERGGG')
    starts.append('')

    # Get all the possible combinations of beginings and endings
    possible_combinations = []
    for start in starts:
        for end in ends:

            possible_combinations.append((start + end).strip(' '))


    # Try to puzzle the Sequences together from the starts and ends
    combinations = 0
    SNPs = 0
    combinations_SNPs = 0
    for seq in seqs:
        is_SNP = False
        # Check if there is an exact match
        if seq in possible_combinations:
            combinations += 1
        else:
            # check if it is an SNP to one of the top combinations
            for top_seq in top_seqs:
                if lstn.distance(seq,top_seq) < 2:

                    SNPs += 1
                    is_SNP = True
                    break
            # check if this is an SNP to one of the shuffled combinations
            if not is_SNP:
                print(seq)

                for combination in possible_combinations:
                    if lstn.distance(seq, combination) < 2:
                        combinations_SNPs += 1
                        is_SNP = True
                        break
            if not is_SNP:
                pass




    # print(f"Reads with rare variants (< 1%):\t\t{n - top_n}/{n}\t{(n - top_n)/n}\n")
    # if not n - top_n == 0:
    #     print(f"Of those:\n"
    #           f"Explained through combining:\t\t\t{combinations - top_n}/{n - top_n}\t\t{(combinations - top_n) / (n - top_n)}\n"
    #           f"Explained through SNPs:\t\t\t\t\t{SNPs}/{n - top_n}\t\t{(SNPs) / (n - top_n)}\n"
    #           f"Explained through combining and SNPs:\t{combinations_SNPs}/{n - top_n}\t\t{combinations_SNPs/(n - top_n)}\n"
    #           f"Not explained: \t\t\t\t\t\t\t{n - combinations - SNPs - combinations_SNPs}/{n - top_n}\t\t"
    #           f"{(n - combinations - SNPs - combinations_SNPs)/(n - top_n)}")

    row =  os.path.basename(file.name)

    results.loc[row,'Reads'] = n
    results.loc[row, 'RareVariants'] = n - top_n
    results.loc[row, 'TemplateSwitches'] = combinations - top_n
    results.loc[row, 'SNPs'] = SNPs
    results.loc[row, 'TemplateSwitch+SNPs'] = combinations_SNPs
    results.loc[row, 'NotExplained'] = n - combinations - SNPs - combinations_SNPs


# results.to_csv('/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/allele_variation_causes_numbers_cut_ANY.csv')


