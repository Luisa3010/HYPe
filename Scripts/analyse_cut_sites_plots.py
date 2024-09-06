###### Exploratory script to find out where the positions are where sequences
# are pieced together and if there are any patterns to it, creating some plots
# on the way


import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt


df = pd.read_csv('/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/HYP1_HVD_cut_sites.csv')
id_counts = Counter(df['Id'])

# Look at the sequences with possible different cut sites
# remove those that have all same motifs before and after and cut position
df_overlaps = df.drop_duplicates(subset=[col for col in df.columns if col in ['Allele1CutPos', 'Allele1MotifBefore','Allele2MotifBefore','Id']])
# Look at the length of overlap (number of times a combination of id and allele1 and allele2 appears)
combination_df = df.groupby(['Id', 'Allele1', 'Allele2']).size().reset_index(name= 'count')
combination_df['count'] = combination_df['count'] - 1

print('Counts of Combinations')



# Plot the histogram of the lengths of overlaps
plt.hist(combination_df['count'], bins=range(0, combination_df['count'].max() + 1), edgecolor= 'black')
plt.xlabel('Number of Overlapping Motifs')
plt.ylabel('Frequency')
plt.xticks(range(0, combination_df['count'].max()))

plt.show()


# Plot the histogram of the potential positions of the cuts
plt.hist(df_overlaps['Allele1CutPos'], bins=range(0, df_overlaps['Allele1CutPos'].max() + 1), edgecolor= 'black')
plt.xlabel('Position of cut')
plt.ylabel('Frequency')
plt.xticks(range(0, df_overlaps['Allele1CutPos'].max()))

plt.show()

# Plot the histogram of the potential positions of the cuts
plt.hist(df_overlaps['Allele2CutPos'], bins=range(0, df_overlaps['Allele2CutPos'].max() + 1), edgecolor= 'black')
plt.xlabel('Position of cut')
plt.ylabel('Frequency')
plt.xticks(range(0, df_overlaps['Allele2CutPos'].max()))

plt.show()




# Look at the Sequences with only one possible cut site
single_ids = [item for item, count in id_counts.items() if count == 1]

filtered_df = df[df['Id'].isin(single_ids)]

a1_after_counts = Counter(filtered_df['Allele1MotifAfter'])
a2_after_counts = Counter(filtered_df['Allele2MotifAfter'])

a1_before_counts = Counter(filtered_df['Allele1MotifBefore'])
a2_before_counts = Counter(filtered_df['Allele2MotifBefore'])

print("Motif after the cut")
print(a1_after_counts)
print(a2_after_counts)

print("Motif before the cut")
print(a1_before_counts)
print(a2_before_counts)

print(filtered_df[['Allele1MotifAfter', 'Allele2MotifAfter']].to_string())