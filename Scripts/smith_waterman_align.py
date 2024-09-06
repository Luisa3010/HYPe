###### Helper script to explore how to create the grep_motifs_sw.py script



from minineedle import needle, smith, core
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import argrelextrema, find_peaks
import time
import sys


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
        alignment: smith.SmithWaterman[str] = smith.SmithWaterman(seq, query)
        alignment.align()

        # Get the score and the alignment matrix
        scores.append(alignment.get_score())
        alignment_mats.append(alignment.get_almatrix())

    best_id = np.argmax(scores)
    best_seq = sequences[best_id]
    best_mat = np.matrix(alignment_mats[best_id])
    # get the indices in the matrix
    pos = np.where(best_mat == np.max(scores))[1]

    return best_seq, pos[0]


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


# all observed nucleotide variations of the different motifs
motifs= {
    "TATGAGCGCGGAGGCGGA": "YERGGG",
    "TATGAGCGCGGAGGCGGG": "YERGGG",
    "TATGAACGCGGAGGCGGA": "YERGGG",
    "TATGAGCGTGGAGGCGGA": "YERGGG",
    "TATGAACGCGGAGGCGGG": "YERGGG",
    "AGTAACCGCGGAGGCGGA": "SNRGGG",
    "AGTAACCGCGGAGGCGGA": "SNRGGG",
    "AGTAACCGCGGGGGCGGA": "SNRGGG",
    "AGTAACCGCGGAGGCGGG": "SNRGGG",
    "AGTAACCGCGGGGGCGGG": "SNRGGG",
    "AGTAACCGTGGAGGCGGA": "SNRGGG",
    "AGTGACCGCGGAGAC": "SDRGD",
    "AGTGACCGCGGAGAT": "SDRGD",
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
    "GGTAACCGCGGGGGCGGA": "GNRGGG"
}

motifs_na_list = list(motifs.keys())
# dna sequences that mark the beginning and end of the HVD
start = "GAAAGTGGTAAAAGACCCGGGAGC"
end = "CATAAACACGGAGGTTATGACGAG"





file_path = "/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/ONT_HYP1_TEST.fasta"

with open(file_path, 'r') as file:
    current_seq_id = None
    sequence = ''
    start_time = time.time()

    for line in file:
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Elapsed time: {elapsed_time} seconds")
        start_time = time.time()
        line = line.strip()
        if line.startswith('>'):  # New sequence header
            sequence = ''
            current_seq_id = line.strip()
        else:
            sequence += line.strip()  # Append sequence without newlines

    # determine if the motifs are on this strand or the complementary strand
            seq_vars = []
            #seq_vars.append("AAACCTTACCAAACTCACAAAATTTTCTAATCTTAAACCAAGTTAAAAGTGTCAAAAAACGTTTTAAGAAGCAATGAACGGCAACAATTTGATGCTGTACATTCTTCTGGCGGGCTTTTGCCTGTTCATCTACGGCACTGAAGCCGGAGGATGCAAGCCGGGACCGCGGGGCCGTCCGGTTCACCGGGCCCTGCGGGGAAGATGGGACCGCCCGGCAAATGCGAGAAACCGCCGCCCAACTATGAGAAGTCGCGCCCCCTGGATCCGAGCACCGCAAAAGATAGGCGGCGCCCGCGGCGTTTTCGGTAAAATTTTTTGAATTTATTTTGTTCGGAATGCTTAAAACAAAATAAAAGGATGCCATGAGCGTTGTCAGAGTTGCCCGTGGGGAATATGACAACAAGTGTCCTGCCGGGCCGCCCGGCGATGTCGGCCCTCCCGGACCGCCGGGTCCCAGCGGAGAGCCGGCAAAATGTTCTCCGCCTGAAGGCTATGAAAGTCGTAAAAGACCCGGGAGCTATGAGCGCGGAGGCGGAAGTGACCGCGGAGGCGGACGTGACAATAAGCGCGGAAGTGACCGCGGAGACTATGAGCGCGGAGGCGGAAGTGACCGCGGAGGCGGAAGTAACCGCGGGGGCGGAAGTGACCGCGGAGACTATGACGCGGAGGCGGAAGTAACCGCGGGGGCGGAAGTAACCGCGGAGGCGGAAGTGACCGCGGAGACTATGAGCGCGGAGGCGGAAGTAACCGCGGGGGCGGAAGTGACCGCGGAGACTATGAGCGCGGAGGCGGAAGTAACCGCGGGGGCGGAAGTAACCGCGGGGGGCGGAAGTAACCGCGGAGGCGGACGTGAAGGCGGAGACTATGAGCGCGGAGGCGGAAGTAACCGCGGAGGCGGGCGTGACCGCGGAGACTATGAGCGCGGAGGCGGGCAAAAGCACGGAGGTTATGACGAGCATCATCTTCGCGCTGTGCGCGGGTTCGTTTGAACTTAAATGGAAGGGGGGGGCAGAAATTTTTGATTTTTTTATTTCCAGTTACGAAGATGAACTGGACCTCCCGGACCTCCGGCGGAAGGACCTCCAGCCCCGGATGGCCCCCCGGGGGCTTGTGAATGCAAATATTAAATGATCAGACGTTGAAAGATGGTCGACTGGACTGCCTCTGCCATTGTAA")
            seq_vars.append(sequence)
            seq_vars.append(complement_dna(seq_vars[0]))
            seq_vars.append(seq_vars[0][::-1])
            seq_vars.append(seq_vars[1][::-1])


            # find the start and end of the HVD as well as the direction in which the DNA should be read
            seq1, start_pos = find_in_seqs(seq_vars, start)
            seq2, end_pos = find_in_seqs(seq_vars, end)
            final_aa_seqs = []
            final_na_seqs = []

            if not seq1 == seq2:
                print("Something went wrong, start and end are not on the same sequence. Ommiting sequence")
            else:
                seq = seq1[start_pos : end_pos - len(end) ]

                value_matrix = np.zeros((len(motifs.keys()),len(seq)+1))
                i = 0
                for motif, motif_names in motifs.items():

                    # Make the alignment
                    alignment: smith.SmithWaterman[str] = smith.SmithWaterman(seq, motif)
                    alignment.align()

                    # Get the score and the alignment matrix
                    score = alignment.get_score()
                    alignment_matrix = np.matrix(alignment.get_almatrix())

                    # for testing
                    # for align in sorted([alignment], reverse=True):
                    #         print(align)
                    # sns.heatmap(alignment.get_almatrix())
                    # plt.show()

                    # get normalized scores for the probability of this motif at each position
                    value_matrix[i, : ] = np.max(alignment_matrix, axis = 0) / len(motif)
                    i += 1

                # smooth the overlap values and get the maximum across all motifs at each point
                window_size = 5
                smooth_matrix = np.zeros((len(motifs.keys()),len(seq)-1))
                for i in range(len(motifs.keys())):
                    smooth_matrix[i,2:] = moving_average(value_matrix[i], window_size)
                    # plt.plot(smooth_matrix[i])
                    #plt.plot(value_matrix[i])
                max_fun = np.append(np.max(smooth_matrix, axis = 0),0)
                #plt.plot(max_fun)

                # get the peaks of the maximum funciton
                x_peaks = find_peaks(max_fun)[0]
                y_peaks = max_fun[x_peaks]

                # if the peak is to low consider it noise and remove it
                is_larger = []
                for j in range(len(y_peaks)):
                    if y_peaks[j] < 0.6:
                        is_larger.append(False)
                    else:
                        is_larger.append(True)
                y_peaks = y_peaks[is_larger]
                x_peaks = x_peaks[is_larger]
                # if the peaks are too close to each other remove the first one
                is_far = []
                for j in range(len(x_peaks)-1):
                    if x_peaks[j+1] - x_peaks[j] < 5:
                        is_far.append(False)
                    else:
                        is_far.append(True)
                is_far.append(True)
                y_peaks = y_peaks[is_far]
                x_peaks = x_peaks[is_far]


                #plt.scatter(x = x_peaks, y = y_peaks)


                # Add labels and title
                #plt.xlabel('Position')
                #plt.ylabel('Overlap Score')
                #plt.title('Motif Positions')
                #plt.legend(motifs.values())

                # Show the plot
                #plt.show()

                aa_sequence = []
                na_sequence = []
                prev_x_peak = 0
                has_bad_motifs = False
                # get the motifs with the highest values at the peak
                for j, x_peak in enumerate(x_peaks):
                    # get the ids of the motifs with the highest value around the smoothed peak as well as the positions of the unsmoothed maxima
                    offset = x_peak -2
                    ids = get_all_argmax(value_matrix[:,offset:offset+5], offset = offset)
                    # take the last of the surrounding peaks # todo consider to adjust for next peak and prev peak distance
                    unsmoothed_x_peak = max([ids[k][1]for k in range(len(ids))])

                    # if there is an unambiguous highest peak, that also perfectly fits the na sequence take this motif as the correct one
                    if len(ids) == 1 and value_matrix[ids[0][0],ids[0][1]] == 1.0:
                        na_motif = motifs_na_list[ids[0][0]]
                        na_sequence.append(na_motif)
                        aa_sequence.append(translate(na_motif))

                    else:
                        motif_len = unsmoothed_x_peak - prev_x_peak

                        # allow for some extra nucleic acids to adjust for overlap and then test again which motif fits best
                        motif_region = seq[prev_x_peak:prev_x_peak + 25]
                        best_aligns = []
                        best_score = 0
                        for motif in motifs.keys():

                            # Make the alignment TODO: this is slow, make it less so!
                            alignment: smith.SmithWaterman[str] = smith.SmithWaterman(
                                motif_region, motif)
                            alignment.change_matrix(
                                 core.ScoreMatrix(match = 4, miss = -2, gap = -2))
                            alignment.align()
                            # if the sequences differ by more than three positions trash them - for length 15 this is a score below (15-3) * 4 - 3 * 2 = 42
                            # and for length 18 this is (18-3) * 4 - 3 * 2 = 54
                            if len(motif) == 15:
                                tresh = 42
                            else:
                                tresh = 54
                            align_score = alignment.get_score()
                            if  align_score > best_score and align_score > tresh:
                                best_aligns = [alignment]
                                best_score = alignment.get_score()
                            elif align_score == best_score and align_score > tresh:
                                best_aligns.append(alignment)


                        # if there is no motif that fits well append an X for now
                        if best_score < tresh:
                            has_bad_motif = True
                            na_sequence.append("X")
                            aa_sequence.append("X")
                        else:
                            if len(best_aligns) == 1:
                                na_motif = best_aligns[0].seq2
                                na_sequence.append(na_motif)
                                aa_sequence.append(translate(na_motif))
                            else:
                                # if the aligned motifs differ just a little take the one that is more frequent in the dataset for now - go by position because motif list is ordered
                                na_motif = best_aligns[0].seq2
                                na_sequence.append(na_motif)
                                aa_sequence.append(translate(na_motif))


                    prev_x_peak = unsmoothed_x_peak

                if not has_bad_motifs:
                    final_na_seqs.append(na_sequence)
                    final_aa_seqs.append(aa_sequence)


                print(na_sequence)



                # print(np.argmax(smooth_matrix[:,x_peaks],axis = 1))


