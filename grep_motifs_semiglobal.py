import os
import numpy as np
from scipy.signal import find_peaks
import argparse
from tqdm import tqdm




def semiglobal_matrix(query, ref, match=1, mismatch=-1, gap=-2):
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
        scores.append(alignment_mat[-1,:].max())
        alignment_mats.append(alignment_mat)

    best_id = np.argmax(scores)
    best_seq = sequences[best_id]
    best_mat = np.matrix(alignment_mats[best_id])
    # get the indices in the matrix
    pos = np.where(best_mat == np.max(scores))

    return best_seq, pos[1][0], np.max(scores)


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


def main():

    # get the arguments from the command line
    parser = argparse.ArgumentParser(description="Count occurrences of specific strings in text files.")
    parser.add_argument('-i', '--input', required=True, help='File or folder with the sequences in fasta format')
    parser.add_argument('-o1', '--na_output_file', required=False, default = '', help='Output file to write the deduced DNA motifs')
    parser.add_argument('-o2', '--aa_output_file', required=False, default = '', help='Output file to write the deduced amino acid motifs')
    parser.add_argument('-b', '--allow_bad', required = False, default = False, help= 'Default False, ommits sequences with ambigous motifs. If set to True, ambigous motifs will shown but marked as an X')
    parser.add_argument('-n', '--num_seqs', required = False, default = float('inf'), help= 'Number of sequences after which to automatically end the run' )
    parser.add_argument('-d', '--diff', required = False, default = 0, help = 'Number of bases that are allowed to differ in a motif')

    args = parser.parse_args()

    ext = ''
    if not args.num_seqs == float('inf'):
        ext += '_n' + str(args.num_seqs)
    if args.allow_bad:
        ext += '_bTrue'
    if not args.diff == 1:
        ext += '_diff' + str(args.diff)



    # if the input file is a folder run on every fasta file present
    if os.path.isdir(args.input):

        # create folder for outputs
        output_folder = args.input + '_output'
        try:
            os.mkdir(output_folder)
        except:
            pass

        for file in os.listdir(args.input):
            filename = os.fsdecode(file)
            if filename.endswith(".fasta") or filename.endswith(".fa"):
                print(filename)
                infile_path = args.input + '/' + filename
                if not args.na_output_file == '':
                    outfile1_path = output_folder + '/' + os.path.splitext(filename)[0] +  '_' + args.na_output_file
                    outfile2_path = output_folder + '/' + os.path.splitext(filename)[0] + '_' + args.aa_output_file
                else:
                    outfile1_path = output_folder + '/' + os.path.splitext(filename)[0] + '_' + 'dna' + ext
                    outfile2_path = output_folder + '/' + os.path.splitext(filename)[0] + '_' + 'aa' + ext

                allow_bad_motifs = args.allow_bad
                num_seqs = float(args.num_seqs)
                max_diff = float(args.diff)

                parse_motifs(infile_path, outfile1_path, outfile2_path,
                             allow_bad_motifs, num_seqs, max_diff)

    # if the input is not a directory run for the given file
    else:

        infile_path = args.input
        print(infile_path)
        if not args.na_output_file == '':
            outfile1_path = args.na_output_file
            outfile2_path = args.aa_output_file
        else:
            outfile1_path = os.path.splitext(infile_path)[0] + '_' + 'dna' + ext
            outfile2_path = os.path.splitext(infile_path)[0] + '_' + 'aa' + ext
        allow_bad_motifs = args.allow_bad
        num_seqs = float(args.num_seqs)
        max_diff = float(args.diff)
        parse_motifs(infile_path, outfile1_path, outfile2_path,
                     allow_bad_motifs, num_seqs, max_diff)


def parse_motifs(infile_path, outfile1_path, outfile2_path, allow_bad_motifs, num_seqs, max_diff = 1):

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

    # dna sequences that mark the beginning and end of the HVD
    start = "GAAAGTGGTAAAAGACCCGGGAGC"
    end = "CATAAACACGGAGGTTATGACGAG"

    total = 0
    omitted = 0

    # ensure the output files are clear
    with open(outfile1_path, 'w'), open(outfile2_path, 'w'):
        pass

    with open(infile_path, 'r') as infile, open(outfile1_path, 'a') as outfile1, open(outfile2_path, 'a') as outfile2:
        current_seq_id = None
        sequence = ''
        n = 0

        for line in infile:
            line = line.strip()
            if line.startswith('>'):  # New sequence header
                # if the required number of sequences has been parsed, end the run
                n += 1
                if n > num_seqs:
                    break

                sequence = ''
                current_seq_id = line.strip()
            else:
                sequence += line.strip()  # Append sequence without newlines

        # determine if the motifs are on this strand or the complementary strand
                seq_vars = []
                seq_vars.append(sequence)
                seq_vars.append(complement_dna(seq_vars[0]))
                seq_vars.append(seq_vars[0][::-1])
                seq_vars.append(seq_vars[1][::-1])


                # find the start and end of the HVD as well as the direction in which the DNA should be read
                seq1, start_pos, start_score = find_in_seqs(seq_vars, start)
                seq2, end_pos, end_score = find_in_seqs(seq_vars, end)


                has_bad_motif = False

                # check if the HVD is found correctly
                has_good_scores = start_score / len(start) > 0.75 and end_score / len(end) > 0.75
                if not seq1 == seq2 or (end_pos - len(end) - start_pos) < 0 or not has_good_scores:
                    omitted += 1
                    total += 1
                else:
                    try:
                        seq = seq1[start_pos : end_pos - len(end) ]

                        # get the semi global alignments for each motif for the sequence
                        value_matrix = np.zeros((len(motifs.keys()),len(seq)+1))
                        motifs_ordered = []
                        for i, motif in enumerate(motifs.keys()):
                            # Get the alignment matrix
                            alignment_matrix = semiglobal_matrix(motif, seq)

                            # get normalized scores for the probability of this motif at each position
                            value_matrix[i, :] = alignment_matrix[-1,:] / len(motif)
                            # store the order of the motifs to retrieve ids later
                            motifs_ordered.append(motif)

                        # get the maximum function of all the alignment scores of all motifs
                        max_fun = np.append(np.max(value_matrix, axis = 0), 0)

                        # get the peaks of the maximum function - peaks are required to be at least 14bp apart
                        x_peaks = find_peaks(max_fun, distance = 14)[0]
                        y_peaks = max_fun[x_peaks]

                        aa_sequence = []
                        dna_sequence = []
                        # get the motifs with the highest values at the peak
                        for i, x_peak in enumerate(x_peaks):

                            # if there is a deletion, a second local alignment has to be done because else unwanted behavior can happen where bases from the previous motif are reused
                            observed_len = x_peaks[i] - x_peaks[
                                i - 1] if i > 0 else x_peak
                            if observed_len == 14 or observed_len == 17 and max_diff > 0:
                                # make the local alignment for each motif and get their scores
                                scores = np.zeros(len(motifs_ordered))
                                for j, motif in enumerate(motifs_ordered):
                                    alignment_matrix = semiglobal_matrix(motif,seq[x_peaks[i - 1] - 1 if i > 0 else 0: x_peaks[i]])
                                    scores[j] = alignment_matrix[-1, :].max() / len(motif)

                                # get ids of the max scores
                                ids = np.where(scores == scores.max())[0]

                            # if there is no deletion use the semi global alignment to get the best motifs
                            else:
                                ids = np.where(value_matrix[:, x_peak] == max_fun[x_peak])[0]
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

                            # check if there are candidates left and if all the potential motifs have the same aa sequence
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

                        # check if there might be motifs missing at the end
                        length_difference = len(seq) - x_peaks[-1]
                        if length_difference > max_diff * 2 :
                            dna_sequence.append('X')
                            aa_sequence.append('X')
                            has_bad_motif = True


                        # write header to files
                        if not (not allow_bad_motifs and has_bad_motif):
                            # print(current_seq_id)

                            outfile1.write(current_seq_id + "\n")
                            outfile2.write(current_seq_id + "\n")
                            # write na and aa sequence to files
                            outfile1.write(' '.join(dna_sequence) + "\n")
                            outfile2.write(' '.join(aa_sequence) + "\n")

                        else:
                            omitted += 1

                        total += 1

                    except:
                        omitted += 1
                        total += 1



    print(f"Total sequences:    {total}")
    print(f"Omitted sequences:  {omitted}")
    print(f"Outputs saved to {os.path.dirname(outfile1.name)}")



if __name__ == "__main__":
    main()
