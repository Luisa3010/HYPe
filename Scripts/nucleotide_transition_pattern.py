##### for a fasta input file get the motifs and create outputs that can be used
# to create a sankey diagram - outdated, use visualize_sequences_as_sankey.py
# instead


import pandas as pd


def read_fasta(file_path):
    with open(file_path, 'r') as file:
        sequences = {}
        current_seq_id = None
        sequence = ''

        for line in file:
            line = line.strip()
            if line.startswith('>'):  # New sequence header
                if current_seq_id:
                    sequences[current_seq_id] = sequence
                # Extract the protein ID or some identifier from the header
                current_seq_id = line  # Takes the line as ID
                sequence = ''  # Reset sequence
            else:
                sequence += line.strip()  # Append sequence without newlines

        # Add the last sequence to the dictionary
        if current_seq_id:
            sequences[current_seq_id] = sequence

    return sequences



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


def create_transition_matrix():
    nucleotides = ['a','c','g','t']
    combinations =[]
    for n1 in nucleotides:
        for n2 in nucleotides:
            for n3 in nucleotides:
                combinations.append(n1+n2 + n3)

    return pd.DataFrame(0, index = combinations, columns = combinations)


# fuzzy search for motifs
def is_match(motif, kmer, max_errors = 1):
    if not len(motif) == len(kmer):
        return False
    errors = 0
    for i in range(len(motif)):
        if not motif[i] == kmer[i]:
            errors += 1
            if errors > max_errors:
                return False
    return True





def split_sequence_by_motifs(seq, motifs):
    translated_seq = translate(seq)

    split_seq = []
    i = 0
    while(i < len(translated_seq)):
        # for motif in motifs:
        #     kmer1 = translated_seq[i:i+5]
        #     kmer2 = translated_seq[i:i+6]
        #
        #     if is_match(motif,kmer1):
        #         split_seq.append(seq[i*3:(i+5)*3])
        #         i += len(motif)
        #         n = 0
        #         break
        #
        #     elif is_match(motif, kmer2):
        #         split_seq.append(seq[i*3:(i+6)*3])
        #         i += len(motif)
        #         n = 0
        #         break

        kmers = [translated_seq[i:i+j] for j in range(4,min(9, len(translated_seq)-i+1))]
        is_found = False
        for motif in motifs:
            for kmer in kmers:
                if is_match(motif, kmer) and not is_found:
                    split_seq.append(seq[i * 3:(i + len(kmer)) * 3])
                    i += len(kmer)
                    for s in split_seq:
                        print(translate(s))
                    print(translated_seq[i:])
                    is_found = True
                    break
            if is_found:
                break
        if not is_found:
            i += 1
            split_seq[-1] += seq[i*3:(i+1)*3]




    return split_seq


# HYP1
motifs = [
    "YERGGG",
    "SNRGGG",
    "RDRGD",
    "SDRGD",
    "SDRGGG",
    "RDNKRG",
    "REGGD",
    "RDRGGG",
    "SDRGE",
    "GNRGGG",
    "RDRGE",
    "SDRGSG",
    "GDRGD",
    "RDRGSG",
    "YERRGG",
    "HERGGG",
    "SDNKRG",
    "YKRGGG",
    "SNHGGG"
]
path = "/home/luisa/Documents/Work/Uni Cambridge/Data/HYP1_sanger/HYP1_nucleotide_repeat_only.fasta"


##HYP3
#
# motifs = {
#     "KKYG",
#     "KKYD"
#     "KKGN",
#     "SDEKE",
#     "NEKEEE",
#     "KDEKEEE",
#     "NDEKEEE",
#     "SDEKEKE",
#     "NDEKEKE",
#     "SDEKEEE",
#     "SDDKEKE",
#     "SDDKEEE",
#     "NDGKEKE",
#     "SDQKEEE",
#     "KEEKEEE",
#     "SDKKEEE",
#     "SDENEKEE"
# }
#
# path = "/home/luisa/Documents/Work/Uni Cambridge/Data/HYP3_sanger/HYP3_nucleotide_repeat_only.fasta"


sequences = read_fasta(path)
matrix = create_transition_matrix()
all_split = []

for seq_id, seq in sequences.items():
    split_seq = split_sequence_by_motifs(seq, motifs)
    all_split.extend(split_seq)
    all_split.append("")

    for i in range(len(split_seq)-1):
        end_bases = (split_seq[i][-3:]).lower()
        start_bases = (split_seq[i+1][:3]).lower()
        print(f"{start_bases}-{end_bases}")
        matrix.at[end_bases,start_bases] += 1
        # if end_bases == "ga":
        #     print(f"{end_bases}-{start_bases}\n{translate(split_seq[i])}-{translate(split_seq[i+1])}")


matrix = matrix.loc[(matrix!=0).any(axis=1)]
matrix = ((matrix.T).loc[((matrix.T)!=0).any(axis=1)]).T


print(matrix)


# get the transitions of nucleotide motifs
motif_transition_matrix = pd.DataFrame(0, index = list(set(all_split)), columns = list(set(all_split)))

for i in range(len(all_split)-1):
    if not (all_split[i] == "" or all_split[i+1] == ""):

        motif_transition_matrix.at[all_split[i],all_split[i+1]] += 1



motif_transition_matrix.to_csv("tmp_output.txt")



# # get the input for sankey maker
# nucleotide
j = 0
for i in range(len(all_split)-1):
    if not (all_split[i] == "" or all_split[i+1] == ""):
        print(f"{translate(all_split[i])}\\n{all_split[i]}{j}[1]{translate(all_split[i+1])}\\n{all_split[i+1]}{j+1}")
        j +=1
    else:
        j = 0
#
# # amino acid
# j = 0
# for i in range(len(all_split)-1):
#     if not (all_split[i] == "" or all_split[i+1] == ""):
#         print(f"{translate(all_split[i])}{j}[1]{translate(all_split[i+1])}{j+1}")
#         j +=1
#     else:
#         j = 0
#

