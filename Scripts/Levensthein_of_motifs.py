#### Get the Levensthein distance of the HYP1 G.Pallida motifs and draw a network
# from them



import pandas as pd
import numpy as np
from Levenshtein import distance as levenshtein_distance
import networkx as nx
import matplotlib.pyplot as plt
import Levenshtein
import pandas as pd

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
        # "AGTGACCGCGGAGAC": "SDRGD",
        # "AGTGACCGCGGAGAT": "SDRGD",
        # "AGTAGCCGCGGAGAC": "SSRGD",
        # "CGTGACCGCGGAGAC": "RDRGD",
        "AGTGACCGCGGAGGCGGA": "SDRGGG",
        "AGCGACCGCGGAGGCGGA": "SDRGGG",
        "AGTGACCGCGGAGGCGGG": "SDRGGG",
        "CGTGACAATAAGCGCGGA": "RDNKRG",
        # "CGTGAAGGCGGAGAC": "REGGD",
        "CGTGACCGCGGAGGCGGA": "RDRGGG",
        # "AGTGACCGCGGAGAG": "SDRGE",
        "GGTAACCGCGGAGGCGGG": "GNRGGG",
        "GGTAACCGCGGAGGCGGA": "GNRGGG",
        "GGTAACCGCGGGGGCGGA": "GNRGGG",
        "CGTGACGATCAGCGCGGA": "RDDQRG",
        # "AGATAACCGCGGAGGCGGA": "XSNRGGG"

}

strings = [str(key) for key in motifs.keys()]

# Create an empty graph
G = nx.Graph()

# Add nodes
for s in strings:
    G.add_node(s)

# Compute Levenshtein distances and add edges
for i in range(len(strings)):
    for j in range(i + 1, len(strings)):
        s1 = strings[i]
        s2 = strings[j]
        distance = Levenshtein.distance(s1, s2)
        similarity = 1 / (1 + distance)  # Convert distance to similarity
        G.add_edge(s1, s2, weight=similarity)

# Draw the network
pos = nx.spring_layout(G, weight='weight')  # Spring layout with edge weights
edges = G.edges(data=True)

# Draw nodes
nx.draw_networkx_nodes(G, pos, node_size=500)

# Draw edges with varying thickness based on weight
# nx.draw_networkx_edges(G, pos, edgelist=edges,
#                        width=[d['weight']*10 for (u, v, d) in edges])

# Draw labels
nx.draw_networkx_labels(G, pos, font_size=12, font_family="sans-serif")

# Show plot
plt.title("Network Visualization of Strings Based on Levenshtein Distances")
plt.show()