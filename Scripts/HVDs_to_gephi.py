###### For a given file of fasta sequences, calculate the Levensthein distance
# between each pair of reads and convert it to a matrix and store it as csv file
# that can be used as input for gephi

import pandas as pd
import Levenshtein


path = '/home/luisa/Documents/Work/Uni_Cambridge/Data/Rostochiensis/PACBIO/L19/L19_HYP3_HVD.fasta'
seqs = []
ids = []
with open(path, 'r') as file:
    for line in file:
        if not line.startswith('>'):
            seqs.append(line.strip('\n'))
        else:
            ids.append(line.strip('\n'))

# seqs = [str(key) for key in motifs.keys()]
# Initialize the matrix
distance_matrix = pd.DataFrame(index=seqs, columns=seqs)

# Compute Levenshtein distances
for i in range(len(seqs)):
    for j in range(len(seqs)):
        if i != j:
            distance_matrix.iloc[i, j] = Levenshtein.distance(seqs[i], seqs[j]) ** 10
        else:
            distance_matrix.iloc[i, j] = 0  # Distance to itself is 0

# Initialize an empty list for edges
edges = []

# Iterate through the matrix to extract edges
for source in distance_matrix.index:
    for target in distance_matrix.columns:
        distance = distance_matrix.at[source, target]
        if source != target:  # Exclude self-loops
            # Convert distance to similarity (e.g., inverse distance)
            similarity = 1 / (1 + distance) if distance != 0 else 0
            edges.append((source, target, similarity))

# Convert the edges list to a DataFrame
edges_df = pd.DataFrame(edges, columns=['Source', 'Target', 'Weight'])

# Save the edge list to a CSV file
edges_df.to_csv('rost_HVD_gephi_edges.csv', index=False)

print("Edge list has been successfully created and saved")