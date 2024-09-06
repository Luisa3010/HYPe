#### for G. Rostochiensis get the length distribution of HVDs and plot them

import matplotlib.pyplot as plt


lengths = []
path = '/home/luisa/Documents/Work/Uni_Cambridge/Data/Rostochiensis/PACBIO/L19/L19_HYP3_HVD.fasta'

with open(path, 'r') as file:
    for line in file:
        if not line.startswith('>'):
            lengths.append(len(line.strip('\n')))
            line = line.strip('\n')
            print(f"{line},{lengths[-1]}")



        else:
            id = line.strip('\n')

        # if len(line) > 910 and len(line) < 921:
        #     print(id)
        #     print(line.strip('\n'))


plt.hist(lengths,25, edgecolor = 'black')
plt.title('Length distribution of HYP3 HVDs in G.Rostochiensis')
plt.ylabel('Number of Reads')
plt.xlabel('HVD Length')
plt.show()