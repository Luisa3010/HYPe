###### Exploratory script to find out where the positions are where sequences
# are pieced together and if there are any patterns to it


import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter


# Read the data
df = pd.read_csv('/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/HYP1_HVD_cut_sites.csv')
df = df.drop('Unnamed: 0', axis = 1)


print(df.columns)

# # Plot the positions of the cuts
# plt.hist(df['Allele1CutPos'], bins = 25, edgecolor='black')
# plt.title("Potential Positions for Cuts in the first Allele")
# plt.xlabel("Motif Position in Sequence")
# plt.ylabel("Frequency")
# plt.show()

# plt.hist(df['Allele2CutPos'], bins = 25, edgecolor='black')
# plt.title("Potential Positions for Cuts in the Second Allele")
# plt.xlabel("Motif Position in Sequence")
# plt.ylabel("Frequency")
# plt.show()

# Get the number of overlapping motifs
#  Group together by ID, Allele1 and Allele2
# grouped_df = df.groupby(['Allele1', 'Allele2', 'Id'])
# grouped_counts = grouped_df.count()['Allele1CutPos']
# counts = []
# for item, count in grouped_counts.items():
#     counts.append(count - 1)
#
# print(sum(counts)/len(counts))
# plt.hist(counts, bins = max(counts) - min(counts))
# plt.title('Histogram of overlap lengths for rare variants')
# plt.show()

# Check if it could be only switches inside the 4 clusters - looking at the allele before
# Dict of motifs in which cluster they might be judging from sequence similarity
cluster_dict_before = {
        'YERGGG': 2,
        'SNRGGG': 2,
        'SDRGGG': 2,
        'RDRGGG': 2,
        'RDNKRG': 3,
        'RDRGD': 4,
        'SDRGD': 4,
        'SDRGE': 4,
        'REGGD': 4
}

cluster_dict_after = {
        'YERGGG': 1,
        'SNRGGG': 2,
        'SDRGGG': 3,
        'RDRGGG': 4,
        'RDNKRG': 4,
        'RDRGD': 4,
        'SDRGD': 3,
        'SDRGE': 3,
        'REGGD': 4
}

# Create new columns containing which cluster Allele1MotifBefore and Allele2MotifBefore come from
df['Allele1MotifBeforeCluster'] = df['Allele1MotifBefore'].map(cluster_dict_before)
df['Allele2MotifBeforeCluster'] = df['Allele2MotifBefore'].map(cluster_dict_before)

# Remove all rows where they are not from the same cluster
df_before_cluster = df[df['Allele1MotifBeforeCluster'] == df['Allele2MotifBeforeCluster']]
print(df_before_cluster)


df['Allele1MotifAfterCluster'] = df['Allele1MotifAfter'].map(cluster_dict_after)
df['Allele2MotifAfterCluster'] = df['Allele2MotifAfter'].map(cluster_dict_after)

# Remove all rows where they are not from the same cluster
df_after_cluster = df[df['Allele1MotifAfterCluster'] == df['Allele2MotifAfterCluster']]
print(df_after_cluster)


# print(len(counts))
# get the amount of rows that overlap between the two
overlap = pd.merge(df_after_cluster, df_before_cluster, how='inner')
print(overlap)



# Find unique Ids (that appear only once)
unique_ids = df['Id'].value_counts() == 1

# Filter the DataFrame to keep only rows with unique Ids
unique_df = df[df['Id'].isin(unique_ids[unique_ids].index)]

print(unique_df)
# print(df[['Allele1MotifBefore', 'Allele2MotifBefore']])

# unique_combinations = df.groupby(['Id', 'Allele1', 'Allele2']).filter(lambda x: len(x) == 1)
# print(unique_combinations)

# Get a matrix of the different occuring combinations

matrix = pd.crosstab(unique_df['Allele1MotifAfter'], unique_df['Allele2MotifAfter'])
print(matrix)
matrix = pd.crosstab(unique_df['Allele1MotifBefore'], unique_df['Allele2MotifBefore'])
print(matrix)
