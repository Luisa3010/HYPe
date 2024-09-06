###### from a file or folder of reads, get those that are rare and get their length
# distribution


from collections import Counter
import matplotlib.pyplot as plt


seqs = []
for i in range(1,11):
    path = f'/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_semiglobal/grouped/heads{i}.phred15.400to1400bp.alignscore2000.fulllength_aa_diff0'
    with open(path, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                seqs.append(line.strip('\n'))

allele_counts = Counter(seqs)
n = len(seqs)

# get the most common alleles for this rep
top_alleles = []
for allele, count in allele_counts.most_common(6):
    if count/n > 0.01:
        top_alleles.append(allele)

# Find all the possible starts and ends
starts = []
ends = []
for allele in top_alleles:
    # cut between any motifs

    motifs = allele.split(' ')
    for i in range(len(motifs)):
        starts.append((' ').join(motifs[:i]) + ' ')
        ends.append((' ').join(motifs[i:]))
ends.append('YERGGG')
starts.append('')

# Get all the possible combinations of starts and endings
possible_combinations = []
for start in starts:
    for end in ends:
        possible_combinations.append((start + end).strip(' '))

# Get all the rare variants
rare_variants = []
for seq in seqs:
    if not seq in top_alleles and seq in possible_combinations:
        rare_variants.append(seq)

# Get the lengths of the rare variants
lengths = []
for allele in rare_variants:
    length = len(allele.split(' '))
    lengths.append(length)

# Plot the lengths
plt.hist(lengths, bins = max(lengths) - min(lengths), edgecolor = 'black')
plt.title("Rare Variants Length Distribution of all Heads reps")
plt.xlabel("Length")
plt.ylabel("Frequency")
plt.show()
