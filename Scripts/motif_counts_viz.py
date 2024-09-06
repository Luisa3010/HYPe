import matplotlib.pyplot as plt
import numpy as np

motifs = {
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

def count_motifs(files):
    # create a dictionary to store the motif counts
    motif_counts = {}
    for motif in motifs:
        motif_counts[motif] = 0

    n= 0
    for path in files:
        with open(path, 'r') as file:
            for i, line in enumerate(file):

                if not line.startswith('>'):

                    line_motifs = line.strip('\n').split(' ')
                    for motif in line_motifs:
                        motif_counts[motif] += 1
        n += i
    return motif_counts, n/2




# get the file names
heads_files = []
tails_files = []
for i in range(1, 11):
    heads_files.append(f"/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_semiglobal/grouped/heads{i}.phred15.400to1400bp.alignscore2000.fulllength_dna_diff0")
    tails_files.append(f"/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_semiglobal/grouped/tails{i}.phred15.400to1400bp.alignscore2000.fulllength_dna_diff0")



heads_counts, num_heads_reads = count_motifs(heads_files)
tails_counts, num_tails_reads = count_motifs(tails_files)


# plot everything


# Values from the dictionaries
keys = list(heads_counts.keys())
values1 = list(heads_counts.values())
values2 = list(tails_counts.values())
normalized_heads_counts = [count/num_heads_reads for count in values1]
normalized_tails_counts = [count/num_tails_reads for count in values2]

# Position of the bars on the x-axis
x = np.arange(len(values1))

# Width of the bars
width = 0.35

fig, ax = plt.subplots()
rects1 = ax.bar(x - width/2, normalized_heads_counts, width, label='heads')
rects2 = ax.bar(x + width/2, normalized_tails_counts, width, label='tails')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_xlabel('Motifs')
ax.set_ylabel('Counts per read')
ax.set_title('Motifs in Heads vs Tails')
ax.set_xticks(x)
ax.set_xticklabels([motifs[key] for key in keys])
plt.xticks(rotation=90)
ax.legend()

# Show the plot
plt.tight_layout()
plt.show()