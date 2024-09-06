##### Helper file to write motif triples with counts into a fasta format


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
      "TATGAGCGCGGAGGCGGA":"YERGGG1",
      "TATGAGCGCGGAGGCGGG":"YERGGG2",
      "TATGAACGCGGAGGCGGA":"YERGGG3",
      "TATGAGCGTGGAGGCGGA":"YERGGG4",
      "TATGAACGCGGAGGCGGG":"YERGGG5",
      "AGTAACCGCGGAGGCGGA":"SNRGGG1",
      "AGTAACCGCGGGGGCGGA":"SNRGGG2",
      "AGTAACCGCGGAGGCGGG":"SNRGGG3",
      "AGTAACCGCGGGGGCGGG":"SNRGGG4",
      "AGTAACCGTGGAGGCGGA":"SNRGGG5",
      "CGTGACCGCGGAGAC":"RDRGD1",
      "AGTGACCGCGGAGAC":"SDRGD1",
      "AGTGACCGCGGAGAT":"SCRGD2",
      "AGTGACCGCGGAGGCGGA":"SDRGGG1",
      "AGCGACCGCGGAGGCGGA":"SDRGGG2",
      "AGTGACCGCGGAGGCGGG":"SDRGGG3",
      "CGTGACAATAAGCGCGGA":"RDNKRG1",
      "CGTGAAGGCGGAGAC":"REGGD1",
      "CGTGACCGCGGAGGCGGA":"RDRGGG1",
      "AGTGACCGCGGAGAG":"SDRGE1"
}




# filepath = "/home/luisa/Documents/Work/Uni_Cambridge/grep_motifs/HYP1_motif_variants.txt"
filepath = "/home/luisa/Documents/Work/Uni_Cambridge/grep_motifs/heads_3motifs_out.txt"

with open(filepath, 'r') as file:
    for line in file:
        input_string, count = (line.split(":"))
        # Resulting output
        output = []

        # Iterate through the string and match keys
        i = 0
        while i < len(input_string):
            for motif in motif_variations:
                if input_string.startswith(motif, i):
                    output.append(motif_variations[motif])
                    i += len(motif) - 1
                    break
            i += 1

        # Concatenate the list into a single string
        result = ''.join(output)

        print('>' + result)
        print(input_string)
        print('Count: ' + count[:-1])



# motifs = []
# with open(filepath, 'r') as file:
#     for line in file:
#         motifs.append(line.strip("\n"))
#
# motifs.pop(-1)
#
# for k1 in motif_variations.keys():
#     for k2 in motif_variations.keys():
#         #for k3 in motif_variations.keys():
#             #print(motifs[i] + motifs[j])
#             #print(f">{translate(motifs[i])+translate(motifs[j])}")
#         print(">" + motif_variations[k1]+ motif_variations[k2])
#         print(k1 + k2 )
#


# input_string = "TATGAGCGCGGAGGCGGATATGAGCGCGGAGGCGGACGTGAAGGCGGAGAC"


