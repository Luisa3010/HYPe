#### create heatmaps for the co-occurence of 2 motifs in the different reps


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import warnings
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap



def read_counts(file_path):
    counts_by_files = {}
    counts = {}
    current_file = ''
    with open(file_path, 'r') as f:
        for line in f:
            if not current_file == '' and line[-8:-1] == ".fastq:":
                counts_by_files[current_file] = counts
            if line[-8:-1] == ".fastq:":
                current_file = (line.strip().split(': ')[0]).split('.')[0]
                counts = {}
            else:
                string, count = line.strip().split(': ')
                counts[string] = int(count)
        counts_by_files[current_file] = counts

    return counts_by_files


def get_motifs(string, all_motifs):
    for motif in all_motifs:
        if string[:len(motif)] == motif:
            motif1 = motif
            break
    for motif in all_motifs:
        if string[-len(motif):] == motif:
            motif2 = motif
            break

    return motif1, motif2


# get the counts and the motifs
infile = "/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT/combs_counts.txt"
motifs_file = "/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT/HYP1_motif_variants.txt"

counts_by_file = read_counts(infile)
#
# heads = {}
# tails = {}
# other = {}
#
# for file, counts in counts_by_file.items():
#     if file[:5] == 'tails':
#         tails[file] = counts
#     elif file[:5] == 'heads':
#         heads[file] = counts
#     else:
#         other[file] = counts


# treshold which motifs to show e.g. 1000 means motifs that occure 1/1000 reads are shown as opaque
freq = 100

with open(motifs_file) as file:
    motifs = [line.strip() for line in file]


fig, axes = plt.subplots(2, 5, figsize = (20, 8))

ratio_tables = []
total_tables = []
for i in range(1,11):
    A = np.zeros((len(motifs), len(motifs)))
    B = np.zeros((len(motifs), len(motifs)))
    ratio_table = pd.DataFrame(A, index = motifs, columns = motifs)
    total_table = pd.DataFrame(B, index = motifs, columns = motifs)



    heads = counts_by_file[f'heads{i}']
    tails = counts_by_file[f'tails{i}']
    norm = heads['@'] / tails['@']

    for key in heads.keys():
        # handle division by zero
        try:
            # get the change between the heads and the tails file
            ratio = (heads[key] - tails[key] * norm ) / (heads[key] + tails[key] * norm )
            totals = min(freq*(heads[key] + tails[key]) / (heads['@'] + tails['@']),1)
        except:
            ratio = np.nan
            totals = 0

        motif1, motif2 = get_motifs(key, motifs)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ratio_table[motif2].loc[motif1] = ratio
            total_table[motif2].loc[motif1] = totals


    # # Define custom colormap
    custom_cmap = LinearSegmentedColormap.from_list(
            "custom_cmap", ["#2dc937", "#e7b416", "#cc3232"]
    )

    # heatmap = sns.heatmap(ratio_table, cmap = custom_cmap)
    # plt.tight_layout()
    # plt.show()

    ax = axes[int((i-1)/5),(i-1)%5]
    heatmap = sns.heatmap(ratio_table, cmap = custom_cmap, ax = ax, yticklabels = ((i-1)%5) == 0, xticklabels = i-1 > 4)


    # Overlay with alpha values to show only frequently expressed motif combinations
    for j in range(ratio_table.shape[0]):
        for k in range(ratio_table.shape[1]):
            heatmap.add_patch(plt.Rectangle((k, j), 1, 1, color = 'white',
                                            alpha = 1 - total_table.iloc[j, k]))

    # Add colorbar
    ax.set_title(f'heads{i} vs tails{i}')

    # store the tables for further analysis
    ratio_tables.append(ratio_table.copy())
    total_tables.append(total_table.copy())


plt.suptitle(f'Motif combinations occurring in more than {1/freq *100}% of reads')
plt.tight_layout()
plt.show()

ratio_tables_df = np.array([df.values for df in ratio_tables])
total_tables_df = np.array([df.values for df in total_tables])


# TODO make this work
# Initialize a result DataFrame
result_df = pd.DataFrame(index=ratio_tables_df[0].index, columns=ratio_tables_df[0].columns)

# Check all positions
for i in range(ratio_tables_df.shape[1]):  # Rows
    for j in range(ratio_tables_df.shape[2]):  # Columns
        result_df.iloc[i, j] = np.all(ratio_tables_df[:, i, j] > 0)

print(result_df)
