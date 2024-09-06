##### Exploratory scrippt to try out an aligner to align motifs to one allele

import pandas as pd
from Bio import Align
import seaborn as sns
import matplotlib.pyplot as plt


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





def split_sequence_by_motifs(trans_seq, motifs):
        translated_seq = translate(seq)

        split_seq = []
        i = 0
        while (i < len(translated_seq)):

            kmers = [translated_seq[i:i + j] for j in
                     range(4, min(9, len(translated_seq) - i + 1))]
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
                split_seq[-1] += seq[i * 3:(i + 1) * 3]

        return split_seq



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




# load sequences
seq = "TATGAGCGCGGAGGCGGAAGTGACCGCGGAGGCGGACGTGACCGCGGAGACTATGAACGCGGAGGCGGAAGTAACCGCGGAGGCGGACGTGACCGCGGAGAGTATGAGCGCGGAGGCGGAAGTAACCGCGGAGGCGGAAGTAACCGCGGAGGCGGACGTGACCGCGGAGAGTATGAGCGCGGAGGCGGAAGTAACCGCGGGGGCGGAAGTGACCGCGGAGAGTATGAGCGCGGAGGCGGAAGTAACCGCGGAGGCGGACGTGAAGGCGGAGACTATGAGCGCGGAGGCGGGAGTAACCGCGGGGGCGGAAGTGACCGCGGAGAGTATGAGCGCGGAGGCGGAAGTAACCGCGGAGGCGGAAGTAACCGCGGAGGCGGACGTGACCGCGGAGACTATGAGCGCGGAGGCGGG"
path = "/home/luisa/Documents/Work/Uni Cambridge/Data/HYP1_sanger/HYP1_nucleotide_repeat_only.fasta"
sequences = read_fasta(path)


# split into motifs
seq_by_motifs = split_sequence_by_motifs(seq, motifs)
ids = [seq_id for seq_id, seq in  sequences.items()]

# use distance measure to get all sequences that are closer than a certain treshold
aligner = Align.PairwiseAligner()
aligner.match_score = 1.0  # Score for matching characters
aligner.mismatch_score = -1.0  # Score for mismatching characters
aligner.open_gap_score = -1/6  # Small gap opening penalty
aligner.extend_gap_score = 0  # No gap extension penalty


matrix = pd.DataFrame(0, index = ids, columns = ids)

for seq_id1, seq1 in sequences.items():
    for seq_id2, seq2 in sequences.items():

        norm_factor = 1 / min(len(seq1),len(seq2)) # normalize so shorter sequences dont get penalized
        score = aligner.score(seq1, seq2) * norm_factor
        matrix.at[seq_id1,seq_id2] = score

print(matrix)
sns.heatmap(matrix)
plt.show()


