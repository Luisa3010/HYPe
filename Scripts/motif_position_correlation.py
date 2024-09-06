###### create a matrix with the correlation of the positional distribution for the
# motifs in G. Pallida HYP1


import pandas as pd

motifs = {
        "TATGAGCGCGGAGGCGGA": "YERGGG",
        "TATGAGCGCGGAGGCGGG": "YERGGG",
        "TATGAACGCGGAGGCGGA": "YERGGG",
        "TATGAGCGTGGAGGCGGA": "YERGGG",
        "TATGAACGCGGAGGCGGG": "YERGGG",
        "TATCAGCGCGGAGGCGGA": "YQRGGG",
        "TATCAGCGCGGAGGCGGG": "YQRGGG",
        "AGTAACCGCGGAGGCGGA": "SNRGGG",
        "AGTAACCGCGGGGGCGGA": "SNRGGG",
        "AGTAACCGCGGAGGCGGG": "SNRGGG",
        "AGTAACCGCGGGGGCGGG": "SNRGGG",
        "AGTAACCGTGGAGGCGGA": "SNRGGG",
        "AGTGACCGCGGAGAC": "SDRGD",
        "AGTGACCGCGGAGAT": "SDRGD",
        "AGTAGCCGCGGAGAC": "SSRGD",
        "CGTGACCGCGGAGAC": "RDRGD",
        "AGTGACCGCGGAGGCGGA": "SDRGGG",
        "AGCGACCGCGGAGGCGGA": "SDRGGG",
        "AGTGACCGCGGAGGCGGG": "SDRGGG",
        "CGTGACAATAAGCGCGGA": "RDNKRG",
        "CGTGAAGGCGGAGAC": "REGGD",
        "CGTGACCGCGGAGGCGGA": "RDRGGG",
        "AGTGACCGCGGAGAG": "SDRGE",
        "GGTAACCGCGGAGGCGGG": "GNRGGG",
        "GGTAACCGCGGAGGCGGA": "GNRGGG",
        "GGTAACCGCGGGGGCGGA": "GNRGGG",
        "CGTGACGATCAGCGCGGA": "RDDQRG",
        # "AGATAACCGCGGAGGCGGA": "XSNRGGG"

}

counts = {}
for motif in motifs.keys():
    counts[motif] = [0 for i in range(40)]

for k in range(1,11):
    path = f"/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_SW/HYP1_ONT_fasta/tails{k}.phred15.400to1400bp.alignscore2000.fulllength_dna"


    with open(path, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                line = line.strip('\n').split(' ')
                for i,motif in enumerate(line):
                    (counts[motif])[i] += 1



for motif, count_list in counts.items():
    relevant_count_list = count_list.copy()
    for i in range(len(count_list)):
        if count_list[i]/max(count_list) > 0.01 and count_list[i] > 1:
            relevant_count_list[i] = count_list[i]
        else:
            relevant_count_list[i] = 0
    counts[motif] = relevant_count_list



# get correlation matrix between different motifs
df = pd.DataFrame(counts)
correlation_matrix = df.corr()
print(correlation_matrix)
df.corr().to_csv('correlation_matrix.csv')

# get highly correlated motifs
labels = [(i, j) for i, j in correlation_matrix[correlation_matrix > 0.95].stack().index.tolist() if i < j]


print(counts['AGTGACCGCGGAGAC'])
print(counts['AGTAGCCGCGGAGAC'])
print(labels)
