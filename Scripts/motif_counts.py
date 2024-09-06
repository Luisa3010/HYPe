from collections import Counter
files = []

for i in range(1, 11):
    # files.append(f"/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_semiglobal/heads{i}.phred15.400to1400bp.alignscore2000.fulllength_dna_diff0")
    files.append(f"/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_semiglobal/tails{i}.phred15.400to1400bp.alignscore2000.fulllength_dna_diff0")


motifs = []
for file in files:
    with open(file, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                motifs.extend(line.strip('\n').split(' '))

motif_counts = Counter(motifs)

for item, count in motif_counts.items():
    print(f"{item} {count}")