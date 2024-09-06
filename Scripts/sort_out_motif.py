##### Sort the Sequences by if they contain the specified motif or not into one
# of two output files


import os

# Path to directory
path = '/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_semiglobal/grouped'
motif = 'AGTGACCGCGGAGAT'

# Iterate over all files
for i in range(1,11):
    filepath = os.path.join(path, f'heads{i}.phred15.400to1400bp.alignscore2000.fulllength_dna_diff0')
    with (open(filepath, 'r') as file,
          open(filepath + f'_{motif}', 'w') as outfile,
          open(filepath + f'_no{motif}', 'w') as no_outfile):
        for line in file:
            if line.startswith('>'):
                id = line
            else:
                # Sort the Sequences by if they contain the motif or not into one of two output files
                if motif in line:
                    outfile.write(id)
                    outfile.write(line)
                else:
                    no_outfile.write(id)
                    no_outfile.write(line)



