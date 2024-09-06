##### Create a csv file for gephi to create a network from the Levensthein distance
# of motifs

import pandas as pd
import Levenshtein

# Example list of strings
# List of strings
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

strings = [str(key) for key in motifs.keys()]
# Initialize the matrix
distance_matrix = pd.DataFrame(index=strings, columns=strings)

# Compute Levenshtein distances
for i in range(len(strings)):
    for j in range(len(strings)):
        if i != j:
            distance_matrix.iloc[i, j] = Levenshtein.distance(strings[i], strings[j]) ** 4
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
edges_df.to_csv('gephi_edges.csv', index=False)

print("Edge list has been successfully created and saved as 'gephi_edges.csv'")