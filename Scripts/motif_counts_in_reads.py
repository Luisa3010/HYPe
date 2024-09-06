###### Get the aggregated counts of motifs for the given input files

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

# initialize all zero dicts to store the counts of motifs
dna_counts = {}
aa_counts = {}
read_count = 0

for key, value in motifs.items():
    dna_counts[key] = 0
    aa_counts[value] = 0

for i in range(1,10):
    path = f"/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_SW/HYP1_ONT_fasta/heads{i}.phred15.400to1400bp.alignscore2000.fulllength_dna"

    # count the occurences of the motifs on dna level
    with open(path, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                line = line.strip('\n').split(' ')
                for dna_motif in motifs.keys():
                    dna_counts[dna_motif] += line.count(dna_motif)
            else:
                read_count += 1


# aggregate the counts on dna level to get the counts on aa level
for key, value in motifs.items():
    aa_counts[value] += dna_counts[key]

# print the results
for key, value in aa_counts.items():
    print(f'{key}, {value}')

print('\n\n')

for key, value in dna_counts.items():
    print(f'{key}, {value}')

print('\n\n')

print(f"Total reads: {read_count}")
