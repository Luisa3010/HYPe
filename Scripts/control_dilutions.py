##### Create a table for the control dilutions, how much of each of the reference
# sequences was actually found in the different mixes


import pandas as pd
import numpy as np
from tqdm import tqdm
import itertools


ref_seqs = ['YERGGG', 'YERGGG SNRGGG RDRGD YERGGG SNRGGG SDRGE YERGGG SNRGGG SDRGD YERGGG SNRGGG RDRGD YERGGG',
            'YERGGG SDRGGG RDNKRG SDRGD YERGGG SDRGGG SNRGGG SDRGD YERGGG SNRGGG SNRGGG SDRGD YERGGG SNRGGG SDRGD YERGGG SNRGGG SNRGGG SNRGGG REGGD YERGGG SNRGGG RDRGD YERGGG',
            'YERGGG SDRGGG RDRGD YERGGG SNRGGG SNRGGG SNRGGG SNRGGG RDRGD YERGGG SNRGGG SDRGD YERGGG SDRGGG RDRGD YERGGG SNRGGG SNRGGG REGGD YERGGG SNRGGG RDRGD YERGGG'
            ]
reps = [1,2,3]
mixs = [1,2]
dils = [0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 128]


path = '/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_semiglobal/'



# Generate all combinations of mix, rep, and dil
combinations = list(itertools.product(mixs, reps, dils))

# Create DataFrame with columns for each of the lists
df = pd.DataFrame(combinations, columns=['Mix', 'Rep', 'Dil'])

# Add four additional columns initialized with zeros
df['Ref1'] = 0
df['Ref2'] = 0
df['Ref3'] = 0
df['Ref4'] = 0
df['Total'] = 0

for mix in mixs:
    print(f"Starting mix {mix}/{len(mixs)}")

    for rep in reps:
        print(f"Starting Rep {rep}/{len(reps)}")

        for dil in tqdm(dils):
            infile = path + f'plasmid_mix{mix}_rep{rep}_dil{dil}in1000.phred15.400to1400bp.alignscore2000.fulllength_aa_diff0'
            with open(infile, 'r') as file:
                n = 0
                for line in file:
                    if not line.startswith('>'):
                        n += 1
                        for i, ref_seq in enumerate(ref_seqs):
                            line = line.strip('\n')
                            line = line.replace('YQRGGG', 'YERGGG').replace('RDDQRG', 'RDNKRG').replace('GNRGGG', 'SNRGGG').replace('SSRGD','SDRGD')
                            if line == ref_seq:
                                df.loc[
                                    (df['Mix'] == mix) &
                                    (df['Rep'] == rep) &
                                    (df['Dil'] == dil),
                                    f'Ref{i + 1}'
                                ] += 1
                                break
                df.loc[
                    (df['Mix'] == mix) &
                    (df['Rep'] == rep) &
                    (df['Dil'] == dil),
                    'Total'
                ] = n
print(df)

print("Saving results to /home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/controls_seqs.csv")
df.to_csv('/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/controls_seqs.csv')