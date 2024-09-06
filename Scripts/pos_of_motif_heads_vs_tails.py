##### for the heads and tails file in G. Pallida HYP1 get the positions of motifs
# and create plots for each motif comparing the positions


import matplotlib.pyplot as plt
from collections import Counter
from tqdm import tqdm

def get_positions(infile, motif):
    abs_positions = []
    rel_positions = []
    with open(infile, "r") as file:
        n = 0
        for i,line in enumerate(file):
            if not line.startswith(">"):
                motifs = line.strip("\n").split(" ")
                abs_pos = get_indexes(motifs, motif)
                rel_pos = [pos / len(motifs) for pos in abs_pos]
                abs_positions.extend(abs_pos)
                rel_positions.extend(rel_pos)
            n +=1
    return abs_positions,rel_positions, n


def get_indexes(lst, element):
    return [index for index, value in enumerate(lst) if value == element]



motifs= {
    "TATGAGCGCGGAGGCGGA": "YERGGG1",
    "TATGAGCGCGGAGGCGGG": "YERGGG2",
    "TATGAACGCGGAGGCGGA": "YERGGG3",
    "TATGAGCGTGGAGGCGGA": "YERGGG4",
    "TATGAACGCGGAGGCGGG": "YERGGG5",
    "TATCAGCGCGGAGGCGGA": "YQRGGG1",
    "TATCAGCGCGGAGGCGGG": "YQRGGG2",
    "AGTAACCGCGGAGGCGGA": "SNRGGG1",
    "AGTAACCGCGGGGGCGGA": "SNRGGG2",
    "AGTAACCGCGGAGGCGGG": "SNRGGG3",
    "AGTAACCGCGGGGGCGGG": "SNRGGG4",
    "AGTAACCGTGGAGGCGGA": "SNRGGG5",
    "AGTGACCGCGGAGAC": "SDRGD1",
    "AGTGACCGCGGAGAT": "SDRGD2",
    "AGTAGCCGCGGAGAC": "SSRGD1",
    "CGTGACCGCGGAGAC": "RDRGD1",
    "AGTGACCGCGGAGGCGGA": "SDRGGG1",
    "AGCGACCGCGGAGGCGGA": "SDRGGG2",
    "AGTGACCGCGGAGGCGGG": "SDRGGG3",
    "CGTGACAATAAGCGCGGA": "RDNKRG1",
    "CGTGAAGGCGGAGAC": "REGGD1",
    "CGTGACCGCGGAGGCGGA": "RDRGGG1",
    "AGTGACCGCGGAGAG": "SDRGE1",
    "GGTAACCGCGGAGGCGGG": "GNRGGG1",
    "GGTAACCGCGGAGGCGGA": "GNRGGG2",
    "GGTAACCGCGGGGGCGGA": "GNRGGG3",
    "CGTGACGATCAGCGCGGA": "RDDQRG1",
    # "AGATAACCGCGGAGGCGGA": "XSNRGGG"

}

total_heads = 0
total_tails = 0

for motif, motif_aa in tqdm(motifs.items()):

    all_heads_abs_pos = []
    all_heads_rel_pos = []
    all_tails_abs_pos = []
    all_tails_rel_pos = []

    for i in range(1,11):
        # if i in [7, 9, 5, 8, 1]: # only for testing purpose
        heads_file = f"/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_semiglobal/heads{i}.phred15.400to1400bp.alignscore2000.fulllength_dna_diff0"
        tails_file = f"/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_semiglobal/tails{i}.phred15.400to1400bp.alignscore2000.fulllength_dna_diff0"
        # else:
        #     tails_file = f"/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_semiglobal/heads{i}.phred15.400to1400bp.alignscore2000.fulllength_dna_diff0"
        #     heads_file = f"/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_semiglobal/tails{i}.phred15.400to1400bp.alignscore2000.fulllength_dna_diff0"


        heads_abs_positions, heads_rel_positions, heads_n = get_positions(heads_file, motif)
        tails_abs_positions, tails_rel_positions, tails_n = get_positions(tails_file, motif)

        all_heads_abs_pos.extend(heads_abs_positions)
        all_heads_rel_pos.extend(heads_rel_positions)
        all_tails_abs_pos.extend(tails_abs_positions)
        all_tails_rel_pos.extend(tails_rel_positions)

        total_heads += heads_n
        total_tails += tails_n

    ### plot bar plot for absolute positions
    counter1 = Counter(all_heads_abs_pos)
    counter2 = Counter(all_tails_abs_pos)

    # Prepare data for the bar chart
    elements1 = list(counter1.keys())
    counts1 = list(counter1.values())

    elements2 = list(counter2.keys())
    counts2 = list(counter2.values())

    # Ensure both datasets have bars aligned properly
    all_elements = sorted(set(elements1) | set(elements2))
    counts1_aligned = [counter1.get(elem, 0) for elem in all_elements]
    c2 = [counter2.get(elem, 0) for elem in all_elements]
    counts2_aligned = [count * total_heads/total_tails for count in c2]


    index = [x for x in all_elements if x < 25]


    # Create bar chart
    bar_width = 0.35


    plt.figure()
    plt.bar(index, counts1_aligned[:len(index)], bar_width, label='Heads')
    plt.bar([i + bar_width for i in index], counts2_aligned[:len(index)], bar_width, label='Tails')

    plt.xlabel('Position')
    plt.ylabel('Frequency')
    plt.title(f'Absolute Positions of {motif_aa} in Heads vs Tails')
    plt.xticks([i for i in range(0,25)])
    plt.legend()


    # plt.savefig(f"/home/luisa/Documents/Work/Uni_Cambridge/Diagrams/Absolute Positions of {motif_aa} in Heads vs Tails.png")
    plt.savefig(
        f"/home/luisa/Documents/Work/Uni_Cambridge/Diagrams/Absolute Positions of {motif_aa} in Heads vs Tails.png")
    plt.close()

    # ### plot distribution for relative positions
    # plt.figure()
    # plt.hist(all_heads_rel_pos, bins=100, label='Heads',alpha = 0.5, density=False)
    # plt.hist(all_tails_rel_pos, bins=100, label='Tails', alpha = 0.5, density=False)
    #
    # plt.xlabel('Relative Position')
    # plt.ylabel('Frequency')
    # plt.title(f'Relative Positions of {motif_aa} in Heads vs Tails')
    # plt.xlim([0,1])
    # plt.legend()
    #
    # plt.savefig(f"/home/luisa/Documents/Work/Uni_Cambridge/Diagrams/Relative Positions of {motif_aa} in Heads vs Tails.png")
    # plt.close()