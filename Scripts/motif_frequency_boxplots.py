import os
import matplotlib.pyplot as plt

def count_motifs(file, motifs_dict):
    motif_counts = {}
    for motif in motifs_dict.keys():
        motif_counts[motif] = 0
    total = 0
    for line in file:
        if not line.startswith('>'):
            total += 1
            line = line.strip('\n')
            motifs = line.split(' ')
            for motif in motifs:
                motif_counts[motif] += 1

    return motif_counts, total

def div(list1,list2):
    res = []
    for x, y in zip(list1,list2):
        res.append(x/y)
    return res


# All currently known motifs
motifs= {
    "TATGAGCGCGGAGGCGGA": "YERGGG1",
    "TATGAGCGCGGAGGCGGG": "YERGGG2",
    "TATGAACGCGGAGGCGGA": "YERGGG3",
    "TATGAGCGTGGAGGCGGA": "YERGGG4",
    "TATGAACGCGGAGGCGGG": "YERGGG5",
    "TATCAGCGCGGAGGCGGA": "YQRGGG1",
    "TATCAGCGCGGAGGCGGG": "YQRGGG2",
    "AGTAACCGCGGAGGCGGA": "SNRGGG1",
    "AGTAACCGCGGGGGCGGA": "SNRGGG2",
    "AGTAACCGCGGAGGCGGG": "SNRGGG3",
    "AGTAACCGCGGGGGCGGG": "SNRGGG4",
    "AGTAACCGTGGAGGCGGA": "SNRGGG5",
    "AGTGACCGCGGAGAC": "SDRGD1",
    "AGTGACCGCGGAGAT": "SDRGD2",
    "AGTAGCCGCGGAGAC": "SSRGD1",
    "CGTGACCGCGGAGAC": "RDRGD1",
    "AGTGACCGCGGAGGCGGA": "SDRGGG1",
    "AGCGACCGCGGAGGCGGA": "SDRGGG2",
    "AGTGACCGCGGAGGCGGG": "SDRGGG3",
    "CGTGACAATAAGCGCGGA": "RDNKRG1",
    "CGTGAAGGCGGAGAC": "REGGD1",
    "CGTGACCGCGGAGGCGGA": "RDRGGG1",
    "AGTGACCGCGGAGAG": "SDRGE1",
    "GGTAACCGCGGAGGCGGG": "GNRGGG1",
    "GGTAACCGCGGAGGCGGA": "GNRGGG2",
    "GGTAACCGCGGGGGCGGA": "GNRGGG3",
    "CGTGACGATCAGCGCGGA": "RDDQRG1",
    # "AGATAACCGCGGAGGCGGA": "XSNRGGG"

}


# Initialize dictionaries to strore the counts
heads_counts_per_file = {}
tails_counts_per_file = {}
for motif in motifs.keys():
    heads_counts_per_file[motif] = []
    tails_counts_per_file[motif] = []

heads_totals = []
tails_totals = []

# Iterate over each file
path = '/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_semiglobal/grouped'
for i in range(1,11):
    h_file = f'heads{i}.phred15.400to1400bp.alignscore2000.fulllength_dna_diff0'
    t_file = f'tails{i}.phred15.400to1400bp.alignscore2000.fulllength_dna_diff0'

    with open(os.path.join(path, h_file), 'r') as heads_file, open(os.path.join(path,t_file), 'r') as tails_file:
        # Get the counts for each motif in this file
        heads_counts, h_total = count_motifs(heads_file, motifs)
        tails_counts, t_total = count_motifs(tails_file, motifs)

        # Add the counts to the dicitonaries
        for motif in motifs.keys():
            heads_counts_per_file[motif].append(heads_counts[motif])
            tails_counts_per_file[motif].append(tails_counts[motif])

        # Add the number of lines to the lists
        heads_totals.append(h_total)
        tails_totals.append(t_total)

# Print the results
print(heads_counts_per_file)
print(tails_counts_per_file)


# Make boxplots for each motif
fig, ax = plt.subplots(nrows = 5, ncols = 6, figsize = (15,15), layout="constrained")
for i,motif in enumerate(motifs.keys()):
    ax[int(i/6),i%6].boxplot([div(heads_counts_per_file[motif],heads_totals),
                 div(tails_counts_per_file[motif],tails_totals)
                 ])
    ax[int(i/6),i%6].title.set_text(motifs[motif])
    ax[int(i/6),i%6].set_xticklabels(['Heads','Tails'])

plt.suptitle('Counts of Motifs per File in Heads vs Tails')
plt.show()



