##### This script takes an input file of parsed reads with bad motifs marked with 'X'
# It outputs a plot of where in the HVDs the most errors occur and which motifs are
# most frequent in these positions as well as a convidence intervall how the errors
# should be if the distribution of errors was random

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binom

# path = "/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/heads1_1000seqs_bad_motifs_stripped"
path = "/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/new_motif_reads_for_analysis"

good_seqs_motif_at_pos = [{} for i in range(35)]
bad_seqs_motif_at_pos = [{} for i in range(35)]
total = [0 for i in range(35)]

with open(path, 'r') as file:

    for line in file:
        if line.startswith('>'): continue
        motifs = line.split('\t')


        if not 'X' in motifs:
            for i in range(len(motifs)):
                total[i] += 1

                motifs[i] = motifs[i].strip()
                if not motifs[i] in good_seqs_motif_at_pos[i].keys():
                    good_seqs_motif_at_pos[i][motifs[i]] = 1
                else:
                    good_seqs_motif_at_pos[i][motifs[i]] += 1
        else:
            for i in range(len(motifs)):
                total[i] += 1
                motifs[i] = motifs[i].strip()
                if not motifs[i] in bad_seqs_motif_at_pos[i].keys():
                    bad_seqs_motif_at_pos[i][motifs[i]] = 1
                else:
                    bad_seqs_motif_at_pos[i][motifs[i]] += 1



print(bad_seqs_motif_at_pos)

cut_off = 27
good_seqs_motif_at_pos = good_seqs_motif_at_pos[:cut_off]
bad_seqs_motif_at_pos = bad_seqs_motif_at_pos[:cut_off]

# calculate the percentage of bad motifs
N = sum(total[:cut_off])
k = 0
for pos_dic in bad_seqs_motif_at_pos:
    if 'X' in pos_dic.keys():
        k += pos_dic['X']


# plot the distribution of the errors
# get the 95% confidence interval at eac position
n =  total[: cut_off]
p = k/N
lower = []
upper = []
expected = []

for i in range(len(n)):
    ci_lower, ci_upper = binom.interval(0.95, n[i], p)
    lower.append(ci_lower)
    upper.append(ci_upper)
    expected.append(n[i]*p)

values = [entry.get('X', 0) for entry in bad_seqs_motif_at_pos]

x_labels = []
for entry in bad_seqs_motif_at_pos:
    max_key = max(entry, key=entry.get)
    x_labels.append(max_key)

x = np.arange(len(values))
plt.bar(x, values)
plt.xticks(ticks = x -1,labels = x_labels)
plt.xticks(rotation=45)
plt.plot(x,lower, color = 'grey', linestyle = '--')
plt.plot(x, upper, color = 'grey', linestyle = '--')
plt.plot(x, expected, color = 'k')
plt.ylabel('number of erros')
plt.xlabel('Most common motif at position')
plt.tight_layout()
plt.show()




#
# all_keys = set()
# for entry in bad_seqs_motif_at_pos:
#     all_keys.update(entry.keys())
# all_keys = list(all_keys)
#
# # Initialize a list of lists to hold the values for each key
# values = {key: [] for key in all_keys}
#
# # Populate the values lists
# for entry in bad_seqs_motif_at_pos:
#     for key in all_keys:
#         if key in entry:
#             values[key].append(entry[key])
#         else:
#             values[key].append(0)  # Use 0 if the key is not present in the dictionary
#
# # Generate x values based on the length of the data
# spacing_factor = 5  # Factor to introduce spacing
# x = np.arange(len(bad_seqs_motif_at_pos)) * spacing_factor
#
# # Plotting
# fig, ax = plt.subplots(figsize=(28, 6))
#
# # Width of each bar
# bar_width = .5
#
# # Offsets for each key to avoid overlap
# offsets = np.arange(len(all_keys)) * bar_width
#
# # Plot each key as a separate set of bars
# for i, key in enumerate(all_keys):
#     ax.bar(x + offsets[i], values[key], bar_width, label=key)
#
# # Add labels and legend
# ax.set_xlabel('Entry Index')
# ax.set_ylabel('Value')
# ax.set_title('Bar Graph of Dictionary Entries')
# ax.set_xticks(x + bar_width * len(all_keys) / 2)
# ax.set_xticklabels(np.arange(len(bad_seqs_motif_at_pos)))
# ax.legend()
#
# plt.show()
#
#
#
#
#
#
#
#
# all_keys = set()
# for entry in good_seqs_motif_at_pos:
#     all_keys.update(entry.keys())
# all_keys = list(all_keys)
#
# # Initialize a list of lists to hold the values for each key
# values = {key: [] for key in all_keys}
#
# # Populate the values lists
# for entry in good_seqs_motif_at_pos:
#     for key in all_keys:
#         if key in entry:
#             values[key].append(entry[key])
#         else:
#             values[key].append(0)  # Use 0 if the key is not present in the dictionary
#
# # Generate x values based on the length of the data
# spacing_factor = 5  # Factor to introduce spacing
# x = np.arange(len(good_seqs_motif_at_pos)) * spacing_factor
#
# # Plotting
# fig, ax = plt.subplots(figsize=(28, 6))
#
# # Width of each bar
# bar_width = .5
# # Offsets for each key to avoid overlap
# offsets = np.arange(len(all_keys)) * bar_width
#
# # Plot each key as a separate set of bars
# for i, key in enumerate(all_keys):
#     ax.bar(x + offsets[i], values[key], bar_width, label=key)
#
# # Add labels and legend
# ax.set_xlabel('Entry Index')
# ax.set_ylabel('Value')
# ax.set_title('Bar Graph of Dictionary Entries')
# ax.set_xticks(x + bar_width * len(all_keys) / 2)
# ax.set_xticklabels(np.arange(len(good_seqs_motif_at_pos)))
# ax.legend()
#
# plt.show()