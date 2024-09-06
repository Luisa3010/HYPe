######## For G.Pallida HYP3 get counts for motifs and their positions using a
# find and replace approach



import json
from matplotlib import pyplot as plt
from collections import Counter
import numpy as np

# Codon table for protein translation
def translate(seq):
    # include finding the start codon
    table = {
            'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
            'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
            'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
            'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
            'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
            'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
            'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
            'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
            'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
            'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            if codon in table.keys():
                protein += table[codon]
            else:
                return ""
    return protein

motif_variations = {
    "KKYG" : {},
    "KDEKEEE" : {},
    "NDEKEEE" : {},
    "SDEKEKE" : {},
    "NDEKEKE" : {},
    "SDENEKEE" : {},
    "SDEKEEE" : {},
    "SDDKEKE" : {},
    "KKYD" : {},
    "SDEKE" : {},
    "NEKEEE" : {},
    "SDDKEEE" : {},
    "NDGKEKE" : {},
    "SDQKEEE" : {},
    "KEEKEEE" : {},
    "SDKKEEE" : {}
    }


# motif_variations = {
#     "DEKEEE": {},
# 	"KKYGS": {},
# 	"KKYGK": {},
# 	"KKYGN": {},
# 	"DEKEKE": {},
# 	"DENEKEE": {},
# 	"KKYDS": {},
# 	"DDKEKE": {},
# 	"DEKE": {},
# 	"DDKEEE": {},
# 	"EKEEE": {},
# 	"DGKEKE": {},
# 	"DQKEEE": {},
# 	"EEKEEE": {},
# 	"DKKEEE": {}
# }
#
# motif_variations = {
#
#     "KKYGKDEKEEE": {},
#     "KKYGNDEKEEE": {},
# #    "KKYGSDEKEKE": {},
#  #   "KKYGNDEKEKE": {},
#     "KKYGSDENEKEE": {},
#     "KKYGSDEKEEE": {},
#  #   "KKYDSDDKEKE": {},
#  #   "KKYGSDEKE": {},
#     "KKYGNEKEEE": {},
#  #   "KKYGNDGKEKE": {},
#     "KKYGSDQKEEE": {},
#     "KKYGSDDKEEE": {},
#     "KKYGKEEKEEE": {},
#     "KKYGSDKKEEE": {}
# }



# motif_variations = {
#        "KKYG" :{},
#        "KKYG":{}
# }


# motif_variations = {
#         "KKYG": {},
#         "KEEE": {},
#         "SDE": {},
#         "KDE": {},
#         "KEKE": {},
#         "NDE ": {},
#         "KEE": {},
#         "NE": {},
#         "SDD ": {},
#         "KKYD": {},
#         "KE": {},
#         "NDG": {},
#         "SDK": {},
#         "SDQ": {}
# }

#
# motif_variations = {
#         "KEKKYGNDEKE": {},
#         "KEKKYGSDEKE": {},
#         "KEKKYDSDDKE": {},
#         "KEKKYGSDQKE": {},
#         "KEKKYGSDKKE": {},
#         "EEKKYGKDEKE": {},
#         "EEKKYGNDEKE": {},
#         "EEKKYGSDEKE": {},
#         "EEKKYGNDGKE": {},
#         "EEKKYGKEEKE": {},
#         "EEKKYGSDKKE": {}
# }
#
#
# motif_variations = {
# 	"NDEKEEEKKYGK": {},
#     "KDEKEEEKKYGKD": {},
# 	"DEKEEEKKYGN": {},
#     "NEKEEEKKYGK":{},
# 	"KDEKEEEKKYGS": {},
# 	"SDEKEKKYGN": {},
#         "DEKEEEKKYGS": {},
#         "SDEKEKEKKYGS": {},
#         "SDDKEEEKKYGN": {},
#         "NDEKEKEKKYGS": {},
# }
#
#
# motif_variations = {
#         "SDEKE": {},
#         "NEKEEE": {},
#         "SDENEKEE": {},
#
# }


# motif_variations = {
#     "KDEKEEE" : {},
#     "NDEKEEE" : {},
#     "SDEKEKE" : {},
#     "NDEKEKE" : {},
#     "SDEKEEE" : {},
#     "SDDKEKE" : {},
#     "SDDKEEE" : {},
#     "NDGKEKE" : {},
#     "SDQKEEE" : {},
#     "KEEKEEE" : {},
#     "SDKKEEE" : {}
# }

motif_variations = {
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
        "CGTGAAGGCGGAGAC": "REGGD",
        "CGTGACCGCGGAGGCGGA": "RDRGGG",
        "AGTGACCGCGGAGAG": "SDRGE",
        "GGTAACCGCGGAGGCGGG": "GNRGGG",
        "GGTAACCGCGGAGGCGGA": "GNRGGG",
        "GGTAACCGCGGGGGCGGA": "GNRGGG",
        "CGTGACAATAAGCGCGGA": "RDNKRG",
        "CGTGACGATCAGCGCGGA": "RDDQRG",
        # "AGATAACCGCGGAGGCGGA": "XSNRGGG"

}


# seq = "tatgagcgcggaggcggaataaaatatgagcgcggaggcggaaatcttcatatgagcgcggaggcggaaatcttcatcattcattatgagcgcggaggcggaaatcttcatcat"

with open("/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP_all_sanger/HYP_nucleotide_stripped_fasta.txt") as file:
    sequences = [line.rstrip() for line in file]

#
# for seq in sequences:
#
#     # translate the sequence
#     protein = translate(seq.upper())
#
#     # get the position of the start motif and end motifs of the HVD in the translated protein sequence
#     # start sequence of the HVD
#     start_motif = "KKYG"
#     end_motif = "SDENEKEE"
#
#     HVD_start_pos = protein.find(start_motif)
#     HVD_end_pos = 0
#     index = 0
#     while index < len(protein):
#         index = protein.find(end_motif, index)
#         if index == -1:
#             break
#         index += len(end_motif)
#         HVD_end_pos = index
#
#
#     # get only the HVD from the DNA and the Protein
#     protein_HVD = protein[HVD_start_pos: HVD_end_pos]
#     DNA_HVD = seq[HVD_start_pos * 3: HVD_end_pos * 3]
#
#
#     # create a dict that contains the different base variations for each motif, as well as their positions as seen from front and back
#     for motif in motif_variations.keys():
#         index = 0
#         block = 1
#         rev_block = 1
#         while index < len(protein):
#             index = protein_HVD.find(motif, index)
#             # stop if the motif isn't found anymore
#             if index == -1:
#                 break
#             # count the occurences of the start motif to find the block
#             block = protein_HVD[:index + len(motif)].count(start_motif)
#
#             rev_block = protein_HVD[index + len(motif):].count(start_motif)
#             # store the base sequence of the motif in the dict
#             variant = DNA_HVD[index * 3:index * 3 + len(motif) * 3]
#             if not variant in motif_variations[motif].keys():
#                 d = {
#                         "count": 1,
#                         "pos" : [index],
#                         "rev_pos": [len(protein_HVD) - index - len(motif)],
#                         "block": [block],
#                         "rev_block": [rev_block]
#
#                 }
#                 motif_variations[motif][variant] = d
#             else:
#                 motif_variations[motif][variant]["count"] += 1
#                 motif_variations[motif][variant]["pos"].append(index)
#                 motif_variations[motif][variant]["rev_pos"].append(len(protein_HVD) - index - len(motif))
#                 motif_variations[motif][variant]["block"].append(block)
#                 motif_variations[motif][variant]["rev_block"].append(rev_block)
#
#
#             index += len(motif)



for variant, motif in motif_variations.items():
    # i = 1
    # #print(motif)
    # for variant in motif_variations[motif]:
    #     #if len(motif) == 7:
    #
    #         # uncomment this for a neat fasta representation
    if len(motif) == 5:
        print(f">{motif}-{variant} ")
        print(f"{variant.ljust(4*3,'-')}")
        # print('>' + motif)
        # print(variant)
        #print(f"{variant} {motif_variations[motif][variant]['count']}")
        # i += 1



# plot the results
for motif in motif_variations:

    # Create subplots
    fig, ax = plt.subplots(max(len(motif_variations[motif].keys()), 2), 1,
                           figsize = (
                                   8, len(motif_variations[motif].keys()) * 3))
    # Calculate the maximum y-limit
    max_count = 0
    max_len = 0

    counters = []
    for variant in motif_variations[motif]:
        counter = Counter(motif_variations[motif][variant]["pos"])
        counters.append(counter)
        max_count = max(max_count, max(counter.values()))
        max_len = max(max(motif_variations[motif][variant]["pos"]), max_len)


    # Plot bar plots for each variant
    # for i, variant in enumerate(motif_variations[motif]):
    #     counter = counters[i]
    #     xs = list(counter.keys())
    #     ys = list(counter.values())
    #     ax[i].bar(xs, ys)
    #     ax[i].set_title(variant)
    #     ax[i].set_ylim(0, max_count + max_count/20)  # Set the same y-limit for all subplots
    #     ax[i].set_xlim(-max_len/10, max_len + max_len/10)
    #
    #
    # # Set the overall title
    # fig.suptitle(motif, fontsize=16)
    #
    # # Adjust layout to prevent overlap
    # plt.tight_layout(rect=[0, 0, 1, 0.95])
    # fig.subplots_adjust(hspace=0.5)  # Adjust the space between subplots

    # Show the plot
    # plt.show()




#with open("/home/luisa/Documents/Work/Uni Cambridge/Data/HYP3_nucleotide_variations_of_motifs4.json", "w") as outfile:
#    json.dump(motif_variations, outfile)



