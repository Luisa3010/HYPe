####### Look at parts of the conserved region that roughly mimic motifs to get
# rough estimate how many errors to expect in the HVD - result errors are not
# random in the conserved regions, some parts have way more errors than others

import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import Levenshtein as lstn

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

def find_pos(query, ref):
    # get end position with semiglobal alignment
    scores = semiglobal_matrix(query, ref)[-1,:]
    end_pos = np.where(scores == scores.max())[0][0]
    # get start position with semiglobal alignment on reversed sequence
    scores = semiglobal_matrix(query[::-1], ref[::-1])[-1,:]
    start_pos = len(ref) - np.where(scores == scores.max())[0][0]
    return start_pos, end_pos


# conserved1_ref_seq = "GGATGCCATGAGCGTTGTCAGAGTTGCCCGTGGGGAATATGACAACAAGTGTCCTGCCGGGCCGCCCGGCGATGTCGGCCCTCCCGGACCGCCGGGTCCCAGCGGAGAGCCGGCAAAATGTTCTCCGCCTGAAGGCTATGAAAGTCGTAAAAGACCCGGGAGC"
# conserved2_ref_seq = "CATAAACACGGAGGTTATGACGAGCATCATCTTCGCGCTGTGCGCGGGTTCGTTTTGAACTTAAATGGGGAGGGGGCAGAAATTTTTGATTTTTTTATTTCCAGTTACGAAGATGAAGCTGGACCTCCCGGACCTCCGGGCCCGGAAGGACCTCCAGGCCCGGATGGCCCCCCGGGGGC"
conserved1_ref_seq = "TATGACAACAAGTGTCCT" # is from conserved 1 and supposed to look similar to motif
conserved2_ref_seq = "AGATGAAGCTGGACC"# is from conserved 2 and supposed to look similar to motif
# seq = "GGATGCCGTTGGCGTTGTCAGAGTTGCCCGTGGGGAATATGACAACAAGTGTATTGCCGGGCCGCCCGGCGATGTCGGCCCCCCCGGACCGCCGGGTCCCAGCGGAGAGCCGGCGAAATATTATCCGCCTGAAGGCTATGAAAGTGGTAAAAGACCCGGGAGCTATGAGCGCGGAGGCGGAAGTGACCGCGGAGGCGGGCGTGACCGCGGAGACTATGAACGCGGAGGCGGAAGTAACCGCGGAGGCGGAAGTAACCGCGGAGGCGGAAGTAACCGCGGAGGCGGAAGTAACCGCGGAGGCGGACAGTGACCGCGGAGACTATGAGCGCGGAGGCGGAAGTAACCGCGGGGGCGGAAGTGACCGCGGAGACTATGAGCGCGGAGGCGGAAGTGACCGCGGAGGCGGACGTGACCGCGGAGACTATGAGCGCGGAGGCGGAAGTAACCGCGGGGGCGGAAGTAACCGCGGAGGCGGACGTGAAGGCGGAGACTATGAGCGCGGAGGCGGAAGTAACCGCGGAGGCGGACGTGACCGCGGAGACTATGAGCGCGGAGGCGGGCATAAACACGGAGGTTATGACGAGCATCATCTTCGCGCTGTGCGCGGGTTCGTTTTGAAAAATGGGCGGGGGCAGGAATTTTTGATTTTTTATTTTCAGTTACGAAGATGAAGCTGGACCTCCCGGACCTCCGGGCCCTGAAGGACCTCCAGGCCCGGATGGACCGCCGGGGGCTTGCGAATGCAAATATTAAATGATCAGACGTTAAGGTCGACTGACTGCCTCTGCCATTGTAA"

# indices = [0,15,30,45,60,]
# parts = [s[i:j] for i,j in zip(indices, indices[1:]+[None])]

infile = "/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_fasta/plasmid_mix1_rep1_dil0in1000.phred15.400to1400bp.alignscore2000.fulllength.fasta"

c1_mm_proportions = []
c2_mm_proportions = []
with open(infile, 'r') as file:

    for i, line in tqdm(enumerate(file)):
        if not line.startswith('>'):
            seq = line.strip("\n")
            # find region in seq
            start, end = find_pos(conserved1_ref_seq,seq)
            
            # levensthein align to that region to get the score
            c1_mm = lstn.distance(conserved1_ref_seq,seq[start:end])
            c1_mm_proportion = c1_mm / len(conserved1_ref_seq)

            # find region in seq
            start, end = find_pos(conserved2_ref_seq, seq)

            # levensthein align to that region to get the score
            c2_mm = lstn.distance(conserved2_ref_seq, seq[start:end])
            c2_mm_proportion = c2_mm / len(conserved2_ref_seq)

            # append to list of missmatches
            if c1_mm_proportion < 0.25 and c2_mm_proportion < 0.25:
                c1_mm_proportions.append(c1_mm_proportion)
                c2_mm_proportions.append(c2_mm_proportion)


        if i > 0 and i%100 == 0:
            print(f"Ratio of errors in c1\n "
                  f"Mean: {sum(c1_mm_proportions) / len(c1_mm_proportions)}")

            print(f"Ratio of errors in c2\n "
                  f"Mean: {sum(c2_mm_proportions) / len(c2_mm_proportions)}")


            plt.boxplot([c1_mm_proportions, c2_mm_proportions])
            plt.xticks([1,2], ["C1","C2"])
            plt.title(f"Error ratio in conserved region (average of {i} reads)")
            plt.ylabel("Proportion incorrect bases")
            # plt.show()

            plt.savefig("boxplot.png")


print(f"Ratio of errors in c1\n "
      f"Mean: {sum(c1_mm_proportions) / len(c1_mm_proportions)}")

print(f"Ratio of errors in c2\n "
      f"Mean: {sum(c2_mm_proportions) / len(c2_mm_proportions)}")

plt.boxplot([c1_mm_proportions, c2_mm_proportions])
plt.xticks([1, 2], ["C1", "C2"])
plt.title(f"Error ratio in conserved region (average of {i} reads)")
plt.ylabel("% incorrect bases")
# plt.show()

plt.savefig("boxplot.png")