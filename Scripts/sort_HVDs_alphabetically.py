##### sort reads from fasta file alphabetically


import os

directory = '/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_semiglobal/grouped'

for infile in os.listdir(directory):
    if not os.path.isdir(os.path.join(directory,infile)):

        with open(os.path.join(directory,infile), 'r') as file:
            seqs = []
            for line in file:
                if not line.startswith('>'):
                    seqs.append(line)


        outfile = os.path.join(directory, infile + '_sorted')
        with open(outfile, 'w') as out:
            seqs = sorted(seqs, reverse = True)
            out.writelines(seqs)