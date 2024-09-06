###### Plot the length of the HVD for heads vs tails in motifs and blocks


import random

import matplotlib.pyplot as plt
import numpy as np


def get_avg(counts):
    avg = 0
    for i, count in enumerate(counts):
        avg += i * count
    avg /= sum(counts)
    return avg



head_num_motifs = [0 for i in range(1,41)]
head_num_blocks = [0 for i in range(1,21)]
head_total = 0

tail_num_motifs = [0 for i in range(1,41)]
tail_num_blocks = [0 for i in range(1,21)]
tail_total = 0

path = '/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_semiglobal'

# ids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
# sample = random.sample(ids,5)

for i in range(1,11):
    # if i in sample:
    #     infile = path + f"/heads{i}.phred15.400to1400bp.alignscore2000.fulllength_aa_diff0"
    # else:
    #     infile = path + f"/tails{i}.phred15.400to1400bp.alignscore2000.fulllength_aa_diff0"
    infile = path + f"/heads{i}.phred15.400to1400bp.alignscore2000.fulllength_aa_diff0"

    with open(infile, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                num_motifs = line.count(' ') + 1
                num_blocks = line.count('YERGGG') + line.count('YQRGGG') - 1
                head_num_motifs[num_motifs] += 1
                head_num_blocks[num_blocks] += 1
                head_total += 1

    # if i in sample:
    #     infile = path + f"/tails{i}.phred15.400to1400bp.alignscore2000.fulllength_aa_diff0"
    # else:
    #     infile = path + f"/heads{i}.phred15.400to1400bp.alignscore2000.fulllength_aa_diff0"
    infile = path + f"/tails{i}.phred15.400to1400bp.alignscore2000.fulllength_aa_diff0"

    with open(infile, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                num_motifs = line.count(' ') + 1
                num_blocks = line.count('YERGGG') + line.count('YQRGGG') - 1
                tail_num_motifs[num_motifs] += 1
                tail_num_blocks[num_blocks] += 1
                tail_total += 1

    
for j in range(0,len(tail_num_motifs)):
    tail_num_motifs[j] /= tail_total
    head_num_motifs[j] /= head_total

for j in range(0,len(tail_num_blocks)):
    tail_num_blocks[j] /= tail_total
    head_num_blocks[j] /= head_total

# print(head_total)
# print(tail_total)

head_avg_motifs = get_avg(head_num_motifs)
head_avg_blocks = get_avg(head_num_blocks)
tail_avg_motifs = get_avg(tail_num_motifs)
tail_avg_blocks = get_avg(tail_num_blocks)

print('head_avg_motifs: ' + str(head_avg_motifs))
print('head_avg_blocks: ' + str(head_avg_blocks))
print('tail_avg_motifs: ' + str(tail_avg_motifs))
print('tail_avg_blocks: ' + str(tail_avg_blocks))


x = np.arange(len(tail_num_motifs))  # the label locations
width = 0.35  # the width of the bars
fig, ax = plt.subplots(figsize = (15,2))
rects1 = ax.bar(x - width/2, head_num_motifs, width, label='Heads')
rects2 = ax.bar(x + width/2, tail_num_motifs, width, label='Tails')
# ax.vlines(head_avg_motifs,ymin = 0, ymax = 1,color = 'k',linestyles = '-')
# ax.vlines(tail_avg_motifs,ymin = 0, ymax = 1,color = 'k',linestyles = '--')
# ax.legend(['Avg Heads','Avg Tails', 'Heads', 'Tails'])
# ax.legend(['Avg Group1','Avg Group2', 'Group1', 'Group2'])
ax.set_ylabel("Proportion of Reads")
ax.set_xlabel("Number of motifs in Read")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)


# ax.set_title('Distribution of HVD Lengths in Motifs in Heads vs Tails')
# ax.set_title('Distribution of HVD Lengths in Motifs in two randomly sampled groups')
# ax.set_yscale('log')


fig.tight_layout()

plt.show()


x = np.arange(len(tail_num_blocks))  # the label locations
width = 0.35  # the width of the bars
fig, ax = plt.subplots()
ax.vlines(head_avg_blocks,ymin = 0, ymax = 1,color = 'k',linestyles = '-')
ax.vlines(tail_avg_blocks,ymin = 0, ymax = 1,color = 'k',linestyles = '--')
rects1 = ax.bar(x - width/2, head_num_blocks, width, label='Heads')
rects2 = ax.bar(x + width/2, tail_num_blocks, width, label='Tails')
ax.legend(['Avg Heads','Avg Tails', 'Heads', 'Tails'])
# ax.legend(['Avg Group1','Avg Group2', 'Group1', 'Group2'])
ax.set_title('Distribution of HVD Lengths in Blocks in Heads vs Tails')
# ax.set_title('Distribution of HVD Lengths in Blocks in two randomly sampled groups')
# ax.set_yscale('log')


fig.tight_layout()

plt.show()