##### Exploratory script to help with the creation of grep_motifs_semiglobal

import os
import numpy as np
from scipy.signal import argrelextrema, find_peaks
import argparse
from statistics import mean
from tqdm import tqdm
import seaborn as sns
import matplotlib.pylab as plt


def semiglobal_matrix(query, ref, match=1, mismatch=-1, gap=-1):
# Fill in  a matrix for semi-global alignment of two sequences
# (A modification of needleman-wunsch)

    nrows=len(query)+1
    ncols=len(ref)+1
    matrix=np.zeros(shape=[nrows, ncols])

    # Whereas the first row needs to be zeros (to allow gaps in front of the query,
    # ...the first column needs successive gap penalties, representing successive
    # deletions of the beginning of the query
    for row in range(1, nrows):
        matrix[row, 0]=matrix[row-1, 0] + gap

    # Filling in the rest of the matrix
    for row in range(1, nrows):
        for col in range(1, ncols):

            # Standard NW procedure: moving from top or left incurs a gap
            topscore=matrix[row-1,col] + gap
            leftscore=matrix[row,col-1] + gap

            # A diagonal move incurs a match or mismatch, depending on the sequence
            if query[row-1] == ref[col-1]:
                diagscore=matrix[row-1,col-1] + match
            else:
                diagscore=matrix[row-1,col-1] + mismatch

            # Each cell is filled with the best-scoring route
            matrix[row,col]=max(topscore, leftscore, diagscore)

    # In the resulting matrix, the last row (minus the first column) gives the best
    # alignment ending at each position in the ref sequence
    return(matrix)



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

    if seq == "AGATAACCGCGGAGGCGGA": return "XSNRGGG"
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            if codon in table.keys():
                protein += table[codon]
            else:
                return ""
    return protein



def complement_dna(sequence):
    # get complementary DNA strand
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complementary_sequence = ''.join([complement[base] for base in sequence])

    return complementary_sequence



def find_in_seqs(sequences, query):
    # returns the position of a subsequence in a number of DNA sequence as well as the sequence where it i most likely to be
    scores = []
    alignment_mats = []
    # get the best alignment scores for the strand, reversed, and complementary
    for seq in sequences:
        # Make the alignment
        alignment_mat = semiglobal_matrix(query, seq)

        # Get the score and the alignment matrix
        scores.append(alignment_mat.max())
        alignment_mats.append(alignment_mat)

    best_id = np.argmax(scores)
    best_seq = sequences[best_id]
    best_mat = np.matrix(alignment_mats[best_id])
    # get the indices in the matrix
    pos = np.where(best_mat == np.max(scores))

    # ax = sns.heatmap(best_mat)
    # plt.show()

    return best_seq, pos[1][0]


def moving_average(data, window_size):
    # smooth the function
    return np.convolve(data, np.ones(window_size) / window_size, mode='valid')


def get_all_argmax(matrix, offset = 0):
    max_value = np.max(matrix)
    ids = list(zip(*np.where(matrix == max_value)))
    for i in range(len(ids)):
        ids[i] = (ids[i][0], ids[i][1] + offset)
    return ids

def to_string(char_list):
    list_as_string = ""
    for char in char_list:
        list_as_string += char
    return list_as_string

def has_all_equal(list):
    e0 = list[0]
    for i in range(len(list)):
        if not list[i] == e0:
            return False
    return True


    # all observed nucleotide variations of the different motifs
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

motifs_na_list = list(motifs.keys())
# dna sequences that mark the beginning and end of the HVD
start = "GAAAGTGGTAAAAGACCCGGGAGC"
end = "CATAAACACGGAGGTTATGACGAG"


# infile_path = "/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/ONT_HYP1_TEST.fasta"
# outfile1_path = "/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/na_TEST.fasta"
# outfile2_path = "/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/aa_TEST.fasta"
# allow_bad_motifs = False

max_diff = 1
allow_bad_motifs = False

# determine if the motifs are on this strand or the complementary strand
seq_vars = []
# seq_vars.append("AAACCTTACCAAACTCACAAAATTTTCTAATCTTAAACCAAGTTAAAAGTGTCAAAAAACGTTTTAAGAAGCAATGAACGGCAACAATTTGATGCTGTACATTCTTCTGGCGGGCTTTTGCCTGTTCATCTACGGCACTGAAGCCGGAGGATGCAAGCCGGGACCGCGGGGCCGTCCGGTTCACCGGGCCCTGCGGGGAAGATGGGACCGCCCGGCAAATGCGAGAAACCGCCGCCCAACTATGAGAAGTCGCGCCCCCTGGATCCGAGCACCGCAAAAGATAGGCGGCGCCCGCGGCGTTTTCGGTAAAATTTTTTGAATTTATTTTGTTCGGAATGCTTAAAACAAAATAAAAGGATGCCATGAGCGTTGTCAGAGTTGCCCGTGGGGAATATGACAACAAGTGTCCTGCCGGGCCGCCCGGCGATGTCGGCCCTCCCGGACCGCCGGGTCCCAGCGGAGAGCCGGCAAAATGTTCTCCGCCTGAAGGCTATGAAAGTCGTAAAAGACCCGGGAGCTATGAGCGCGGAGGCGGAAGTGACCGCGGAGGCGGACGTGACAATAAGCGCGGAAGTGACCGCGGAGACTATGAGCGCGGAGGCGGAAGTGACCGCGGAGGCGGAAGTAACCGCGGGGGCGGAAGTGACCGCGGAGACTATGACGCGGAGGCGGAAGTAACCGCGGGGGCGGAAGTAACCGCGGAGGCGGAAGTGACCGCGGAGACTATGAGCGCGGAGGCGGAAGTAACCGCGGGGGCGGAAGTGACCGCGGAGACTATGAGCGCGGAGGCGGAAGTAACCGCGGGGGCGGAAGTAACCGCGGGGGGCGGAAGTAACCGCGGAGGCGGACGTGAAGGCGGAGACTATGAGCGCGGAGGCGGAAGTAACCGCGGAGGCGGGCGTGACCGCGGAGACTATGAGCGCGGAGGCGGGCAAAAGCACGGAGGTTATGACGAGCATCATCTTCGCGCTGTGCGCGGGTTCGTTTGAACTTAAATGGAAGGGGGGGGCAGAAATTTTTGATTTTTTTATTTCCAGTTACGAAGATGAACTGGACCTCCCGGACCTCCGGCGGAAGGACCTCCAGCCCCGGATGGCCCCCCGGGGGCTTGTGAATGCAAATATTAAATGATCAGACGTTGAAAGATGGTCGACTGGACTGCCTCTGCCATTGTAA")
seq_vars.append("AAACCTTACCAAACTCACAAAATTTTCTAATCTTAAACCAAGTTAAAAGTGTCAAAAAACGTTTTAAGAAGCAATGAACGGCAACAATTTGATGCTGTACATTCTTCTGGCGGGCTTTTGCCTGTTCATCTACGGCACTGAAGCCGGAGGATGCAAGCCGGGACCACAGGGGCCGTCCGGTTCACCGGGCCCTGCGGGGAAGATGGGACCGCCCGGCAAATGCGAGAAACCGCCGCCCAACTATGAGAAGTCGCGCCCCCTAAATACGAGCACCGCAAAAGATCAGCGGCGCCCGCGGCGTTTTCGGTAAAATTTTTTGAATTTATTTTTTTTTGTTCGGAATGCTTAAAACAAAATAAAAGGATGCCATGAGCGTTGTCAGAGTTGCCCGTGGGGAATATGACAACAAGTGTCCTGCCGGGCCGCCCGGCGATGTCAGGCCCCCCCGGACCGCCGGGTCCCAGCGGAGAGCCGGCAAAATGTTCTCCGCCTGAAGGCTATGAAAGTGGTAAAAGACCCGGGAGCTATGAGCGCGGAGGCGGAAGTAACCGCGGAGGCGGGCGTGACCGCGGAGACTATGAGCGCGGAGGCGGAAGTAACCGCGGGGGCGGAAGTGACCGCGGAGAGTACGAGCGCGGAGGCGGATGACTGCGGAGGCGGAAGTGACCGCGGAGATTATGAGCGCGGAGGCGGAAGTAACCGCGGAGGCGGGCGTGACCGCGGAGACTATGAGCGCGGAGGCGGGCATAAACACGGAGGTTATGACGAGCATCATCTTCGCGCTGTGCGCGGGTTCGTTTGAACTTAAATGGAAGGGGGGGGCAGAAATTTTTGATTTTTTTATTTCCAGTTACGAAGATGAAGCTCGGACCTCCCGGACCTCCGGGCCCGGAAGGACCTCCAGGCCCGGATGGCCCCCGGGGGCTTGTGAATGCAAATATTAAATGATCAGACGTTGAAAGATGGTCGACTGGACTGCCTCTGCCATTGTAA")
# seq_vars.append(sequence)
seq_vars.append(complement_dna(seq_vars[0]))
seq_vars.append(seq_vars[0][::-1])
seq_vars.append(seq_vars[1][::-1])


# find the start and end of the HVD as well as the direction in which the DNA should be read
seq1, start_pos = find_in_seqs(seq_vars, start)
seq2, end_pos = find_in_seqs(seq_vars, end)

seq = seq1[start_pos : end_pos - len(end)]
value_matrix = np.zeros((len(motifs.keys()),len(seq)+1))
i = 0

motifs_ordered = []
for motif, motif_names in motifs.items():
    # Get the alignment matrix
    alignment_matrix = semiglobal_matrix(motif, seq)

    # get normalized scores for the probability of this motif at each position
    value_matrix[i, :] = alignment_matrix[-1,:] / len(motif)
    # ax = sns.heatmap(alignment_matrix)
    # plt.show()
    motifs_ordered.append(motif)
    i += 1

# get the maximum function of all the scores of all motifs
max_fun = np.append(np.max(value_matrix, axis = 0),0)

# get the peaks of the maximum function - peaks are required to be at least 12bp apart
min_dist = 12
x_peaks = find_peaks(max_fun, distance = min_dist)[0]
y_peaks = max_fun[x_peaks]


plt.figure(figsize=(20,10))
plt.plot(value_matrix.T)
plt.plot(max_fun, color = 'black', linewidth = 2)
plt.scatter(x_peaks,y_peaks, color = 'r')
legend = list(motifs.values())
legend.append('Maximum')
plt.legend(legend,bbox_to_anchor=(1.1, 1))
plt.show()



aa_sequence = []
dna_sequence = []
prev_x_peak = 0
# get the motifs with the highest values at the peak

# TODO: add this to script!
for i, x_peak in enumerate(x_peaks):

    # if there is a deletion, a second local alignment has to be done because else unwanted behavior can happen where bases from the previous motif are reused
    observed_len = x_peaks[i] - x_peaks[i-1] if i > 0 else x_peak
    print(observed_len)
    if observed_len == 14 or observed_len == 17 and max_diff > 0:
    # if True:
        # make the local alignment for each motif and get their scores
        scores = np.zeros(len(motifs_ordered))
        for j, motif in enumerate(motifs_ordered):
            alignment_matrix = semiglobal_matrix(motif, seq[x_peaks[i-1]-1 if i > 0 else 0:x_peaks[i]],1,-1,-1)
            scores[j] = alignment_matrix[-1,:].max()/len(motif)

        # get ids of the max scores
        ids = np.where(scores == scores.max())[0]

    # if there is no deletion use the semi global alignment to get the best motifs
    else:
        ids = np.where(value_matrix[:,x_peak] == max_fun[x_peak])[0]
        scores = value_matrix[:, x_peak]

        # check all the motifs with the best aligning sequences at this position
    potential_motifs = []
    for id in ids:
        dna_motif = motifs_ordered[id]
        # check if the alignment fits well enough
        threshold = (len(dna_motif) - max_diff * 2) / len(dna_motif)
        # check for indels between motifs
        dist = abs(observed_len - len(dna_motif))
        # consider the motif if the allignment is good enough and no to large indels
        if scores[id] >= threshold and dist <= max_diff:
            potential_motifs.append(dna_motif)





    # check if there are candidates left and if all the potential motifs have the same aa sequence, take the most common dna motif
    if len(potential_motifs) > 0 and has_all_equal(potential_motifs):
        dna_motif = potential_motifs[0]
        aa_motif = translate(dna_motif)
        dna_sequence.append(dna_motif)
        aa_sequence.append(aa_motif)
    # if the alignment is not good enough or to large indels append 'X' as a bad motif
    else:
        aa_sequence.append('X')
        dna_sequence.append('X')
        has_bad_motif = True
        if not allow_bad_motifs:
            break


print(aa_sequence)
print(dna_sequence)


    #     # check if there is an insertion between motifs and if there is ommit the sequence
    #     motif = motifs_na_list[ids[0][0]]
    #     if abs(motif_len - len(motif)) > 2:
    #         has_bad_motif = True
    #         # print(f"{j} {motif_region}")
    #     else:
    #         na_motif = motif
    #
    #         na_sequence.append(na_motif)
    #         aa_sequence.append(translate(na_motif))
    #
    # else:
    #     # allow for some extra nucleic acids to adjust for overlap and then test again which motif fits best
    #     motif_region = seq[prev_x_peak:prev_x_peak + 20]
    #     best_aligns = []
    #     best_score = 0
    #     for motif in motifs.keys():
    #
    #         # check if the motif is roughly the correct length
    #         if abs(motif_len - len(motif)) > 2:
    #             continue
    #
    #         # Make another local alignment # TODO make this with a levensthein distance
    #         alignment: smith.SmithWaterman[str] = smith.SmithWaterman(
    #                 motif_region, motif)
    #         alignment.change_matrix(
    #                 # core.ScoreMatrix(match = 4, miss = -2, gap = -2))
    #                 core.ScoreMatrix(match = 4, miss = -2, gap = -2))
    #         alignment.align()
    #         # if the motif differ by more than max_diff from the actual sequence ommit the sequence - this is roughly equal to maxdiff bases difference
    #         max_diff = 1
    #         # max_diff = 0
    #         thresh = (len(motif) - max_diff) * 4 - max_diff * 2
    #         align_score = alignment.get_score()
    #         if align_score > best_score and align_score > thresh:
    #             best_aligns = [alignment]
    #             best_score = alignment.get_score()
    #         elif align_score == best_score and align_score > thresh:
    #             best_aligns.append(alignment)
    #
    #     # if there is no motif that fits well append an X for now
    #     if best_score <= thresh:
    #         has_bad_motif = True
    #         # print(f'{j} {motif_region}')
    #         na_sequence.append("X")
    #         aa_sequence.append("X")
    #         if not allow_bad_motifs:
    #             break
    #     else:
    #         if len(best_aligns) == 1:
    #             na_motif = best_aligns[0].seq2
    #             na_sequence.append(na_motif)
    #             aa_sequence.append(translate(na_motif))
    #         else:
    #             is_same_motif = True
    #             for i in range(len(best_aligns) - 1):
    #                 is_same_motif = is_same_motif and translate(
    #                     best_aligns[i].seq2) == translate(
    #                     best_aligns[i + 1].seq2)
    #             # if the best scoring aligned motifs differ just in the nucleic acid and not the protein take the one that is more frequent in the dataset for now - go by position because motif list is ordered
    #             if is_same_motif:
    #                 na_motif = best_aligns[0].seq2
    #                 na_sequence.append(na_motif)
    #                 aa_sequence.append(translate(na_motif))
    #             # if there is a different motif on protein level, ommit the sequence
    #             else:
    #                 has_bad_motif = True
    #                 # print(f"{j} {motif_region}")
    #                 na_sequence.append("X")
    #                 aa_sequence.append("X")
    #                 if not allow_bad_motifs:
    #                     break











