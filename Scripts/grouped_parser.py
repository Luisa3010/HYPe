##### Take a file of parsed G.Pallida HYP1 motifs and replace those that
# are suposed to be grouped together and store the results in a new file

import os

dir = '/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_semiglobal'
out_dir = '/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_semiglobal/grouped'

files = sorted(os.listdir(dir))


for infile in files[1:]:
    with open(os.path.join(dir,infile), 'r') as file, open(os.path.join(out_dir,infile), 'w') as out:
        for line in file:
            line = line.replace('YQRGGG', 'YERGGG').replace('RDDQRG',
                                                               'RDNKRG').replace(
                'GNRGGG', 'SNRGGG').replace('SSRGD','SDRGD')

            out.write(line)