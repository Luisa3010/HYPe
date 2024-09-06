###### get the difference of occurence for motifs or 2 motif combs in heads and tails
# and visualize as boxplots


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

def read_counts(file_path):
    counts = {}
    with open(file_path, 'r') as f:
        for line in f:
            string, count = line.strip().split(': ')
            counts[string] = int(count)
    return counts

def compare_counts(file1, file2, output_image='comparison_chart.png'):
    counts1 = read_counts(file1)
    counts2 = read_counts(file2)

    # Ensure both dictionaries have the same keys
    all_strings = sorted(set(counts1.keys()).union(set(counts2.keys())))
    counts1 = {string: counts1.get(string, 0) for string in all_strings}
    counts2 = {string: counts2.get(string, 0) for string in all_strings}

    norm = counts1["@"]/counts2["@"]
    # Create a DataFrame for easier plotting
    df = pd.DataFrame({
        'String': list(all_strings),
        'Count1': [counts1[string] for string in all_strings],
        'Count2': [counts2[string] for string in all_strings]
    })



    # Plotting the bar chart
    x = range(len(df))
    # plt.figure(figsize=(10, 6))
    # # plt.bar(x, df['Count1'], width=0.4, label='Heads', align='center')
    # # plt.bar(x, df['Count2'] * norm , width=0.4, label='Tails', align='edge')
    # plt.bar(x, df['Count1'] / df['Count2'] / norm) # check if the motif is up or down regulated
    # plt.xticks(ticks=x, labels=df['String'], rotation=90)
    # plt.xlabel('Motif Variant')
    # plt.ylabel('Counts')
    #
    # plt.tight_layout()
    # plt.savefig(output_image)
    # plt.show()
    return df['String'], df['Count2'] / df['Count1'] * norm

# file1 = '/home/luisa/Documents/Work/Uni_Cambridge/grep_motifs/heads_out.txt'
# file2 = '/home/luisa/Documents/Work/Uni_Cambridge/grep_motifs/tails_out.txt'
#
#
# control_file1 = '/home/luisa/Documents/Work/Uni_Cambridge/grep_motifs/control_heads1_out.txt'
# control_file2 = '/home/luisa/Documents/Work/Uni_Cambridge/grep_motifs/control_heads2_out.txt'

dir = '/home/luisa/Documents/Work/Uni_Cambridge/grep_motifs/'


diffs = np.zeros((10,21))
for i in range(10):
    motifs, diff = compare_counts(f"{dir}heads{i+1}_out.txt",f"{dir}tails{i+1}_out.txt")
    diffs[i,:] = diff

plt.figure(figsize=(10, 6))
sns.boxplot(data = diffs[:,1:21])
plt.xticks(ticks=range(20), labels=motifs[1:21], rotation=90)
plt.tight_layout()
plt.show()



diffs = np.zeros((10,401))
for i in range(10):
    motifs, diff = compare_counts(f"{dir}comb_heads{i+1}_out.txt",f"{dir}comb_tails{i+1}_out.txt")
    diffs[i,:] = diff

motifs = motifs[~np.isnan(diffs).all(axis=0)]
diffs = diffs[:, ~np.isnan(diffs).all(axis=0)]

print(diffs)

step_size = 40
for i in range(0,205,step_size):
    plt.figure(figsize=(10, 6))
    plt.axhline(y = 1, color = 'c', linestyle ='--')

    #sns.boxplot(data = np.log(diffs[:,i:i+step_size]))
    sns.boxplot(data = diffs[:, i:i + step_size])
    plt.xticks(ticks=range(min(step_size,len(diffs[:,i:i + step_size][0]))), labels=motifs[i:i+step_size], rotation=90)
    plt.tight_layout()
    plt.ylim((-0.3,5.3))
    plt.show()





