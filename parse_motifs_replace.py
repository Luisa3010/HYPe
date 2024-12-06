# iteratively replace best matching motifs
import itertools
import os
import warnings
from semiglobal import semiglobal_matrix
import re
import numpy as np
from grep_motifs_semiglobal import translate
import argparse
import math
from itertools import chain, combinations
from tqdm import tqdm


motifs = {
        "TATGAGCGCGGAGGCGGA": "YERGGG",
        "TATGAGCGCGGAGGCGGG": "YERGGG",
        "TATGAACGCGGAGGCGGA": "YERGGG",
        "TATGAGCGTGGAGGCGGA": "YERGGG",
        "TATGAACGCGGAGGCGGG": "YERGGG",
        # "TATCAGCGCGGAGGCGGA": "YQRGGG",
        # "TATCAGCGCGGAGGCGGG": "YQRGGG",
        "AGTAACCGCGGAGGCGGA": "SNRGGG",
        "AGTAACCGCGGGGGCGGA": "SNRGGG",
        "AGTAACCGCGGAGGCGGG": "SNRGGG",
        "AGTAACCGCGGGGGCGGG": "SNRGGG",
        "AGTAACCGTGGAGGCGGA": "SNRGGG",
        "AGTGACCGCGGAGAC": "SDRGD",
        "AGTGACCGCGGAGAT": "SDRGD",
        # "AGTAGCCGCGGAGAC": "SSRGD",
        "CGTGACCGCGGAGAC": "RDRGD",
        "AGTGACCGCGGAGGCGGA": "SDRGGG",
        "AGCGACCGCGGAGGCGGA": "SDRGGG",
        "AGTGACCGCGGAGGCGGG": "SDRGGG",
        "CGTGACAATAAGCGCGGA": "RDNKRG",
        "CGTGAAGGCGGAGAC": "REGGD",
        "CGTGACCGCGGAGGCGGA": "RDRGGG",
        "AGTGACCGCGGAGAG": "SDRGE",
        "CGTGACCGCGGAGAG": "RDRGE"
        # "GGTAACCGCGGAGGCGGG": "GNRGGG",
        # "GGTAACCGCGGAGGCGGA": "GNRGGG",
        # "GGTAACCGCGGGGGCGGA": "GNRGGG",
        # "CGTGACGATCAGCGCGGA": "RDDQRG",
}

# quality scores are calculated as highest score^2 / second_highest score
perfect_match_qs = {
    'TATGAGCGCGGAGGCGGA': 1.125,
    'TATGAGCGCGGAGGCGGG': 1.125,
    'TATGAACGCGGAGGCGGA': 1.125,
    'TATGAGCGTGGAGGCGGA': 1.125,
    'TATGAACGCGGAGGCGGG': 1.125,
    'TATCAGCGCGGAGGCGGA': 1.125,
    'TATCAGCGCGGAGGCGGG': 1.125,
    'AGTAACCGCGGAGGCGGA': 1.125,
    'AGTAACCGCGGGGGCGGA': 1.125,
    'AGTAACCGCGGAGGCGGG': 1.125,
    'AGTAACCGCGGGGGCGGG': 1.125,
    'AGTAACCGTGGAGGCGGA': 1.125,
    'AGTGACCGCGGAGAC': 1.1538461538461537,
    'AGTGACCGCGGAGAT': 1.1538461538461537,
    'AGTAGCCGCGGAGAC': 1.3636363636363638,
    'CGTGACCGCGGAGAC': 1.1538461538461537,
    'AGTGACCGCGGAGGCGGA': 1.125,
    'AGCGACCGCGGAGGCGGA': 1.125,
    'AGTGACCGCGGAGGCGGG': 1.125,
    'CGTGACAATAAGCGCGGA': 1.2857142857142856,
    'CGTGAAGGCGGAGAC': 1.3636363636363638,
    'CGTGACCGCGGAGGCGGA': 1.125,
    'AGTGACCGCGGAGAG': 1.1538461538461537,
    'GGTAACCGCGGAGGCGGG': 1.125,
    'GGTAACCGCGGAGGCGGA': 1.125,
    'GGTAACCGCGGGGGCGGA': 1.125,
    'CGTGACGATCAGCGCGGA': 1.2857142857142856,
    'CGTGACCGCGGAGAG': 1.1538461538461537
}


def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join(complement.get(base, base) for base in reversed(seq))




def semiglobal_matrix(query, ref, match=1, mismatch=-1, gap=-2):
    ''' Fill in  a matrix for semi-global alignment of two sequences
    (A modification of needleman-wunsch) '''

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



def find_in_read(dna, read):
    ''' get the most likely position of a sub-sequence in a read and the normalized alignment score '''

    scores = semiglobal_matrix(dna, read)[-1, :]
    score = max(scores)/len(dna)
    pos = np.argmax(scores)

    return pos, score


def finalize_read(read):
    ''' Transform the finished parsed read into a beautiful form '''

    read = re.sub(r'\s+', ' ', read)

    read = read.strip(' ')
    return read.upper()


def replace_exact_matches(read, motifs):
    ''' Delete instances where there are exact matches of the motifs in the read
    and write their amino acid sequence into a new string. Also generate a
    string containing the quality scores for the replaced motifs '''

    parsed_read = ' ' * len(read)
    quality_string = ' ' ' ' * len(read)

    for dna1, aa in motifs.items():

        while not read.find(dna1) == -1:
            is_overlapping = False
            pos = read.find(dna1)
            if pos == -1:
                continue


            for dna2 in motifs.keys():
                if not dna1 == dna2 and not read[pos - len(dna2) +1 : pos + len(dna1) + len(dna2) -1 ].find(dna2) == -1:
                    is_overlapping = True
                    break

            # If two perfectly matching motifs are not overlapping, remove them from the original read and add the translated motifs to the parsed read
            if not is_overlapping:
                read = read[:pos] + ' ' * len(dna1) + read[pos + len(dna1):]
                # parsed_read = parsed_read[:pos] + '_' * len(aa) + aa + '_' * len(aa) + parsed_read[pos + len(dna1):]
                parsed_read = parsed_read[:pos] + ' ' * len(aa) + aa + ' ' * len(aa) + parsed_read[pos + len(dna1):]
                quality_string = quality_string[:pos] + ' ' * len(aa) + str(round(perfect_match_qs[dna1],2)).ljust(len(aa),' ') + ' ' * len(aa) + quality_string[pos + len(dna1):]
            # if two motifs are overlapping, set to lower case, so it is not find in the next iteration
            else:
                read = read[:pos] + read[pos:pos + len(dna1)].lower() + read[pos + len(dna1):]


    return read.upper(), parsed_read, quality_string


def powerset(iterable):
    '''create the powerset for a given set'''

    s = list(iterable)
    return chain.from_iterable(
            combinations(s, r) for r in range(len(s) + 1))

def has_gaps(config: list):

    # Make sure other cubes don't have overlap
    for i in range(len(config)-1):
        dist = config[i + 1].pos1 - config[i].pos2
        if dist > 1 and dist < 15:
            return True

    return False


def get_possible_configurations(cubes):
    ''' With a given list of cubes generate all configurations of those cube
    where there are no two cubes overlapping '''

    all_configs = list(powerset(cubes))

    possible_configs = []
    for config in all_configs:
        has_overlap = False
        # Remove all configurations where cubes are overlapping
        for i,j in list(itertools.product(range(len(config)),range(len(config)))):
            if (config[i]).is_overlapping(config[j]) and not i == j:
                has_overlap = True
                break

        # Remove all configurations where there is a gap between cubes - gaps that could be full motifs are ok (>= 15bp)

        if not has_overlap and not has_gaps(config):
            possible_configs.append(config)

    return possible_configs


def get_best_config(configs):
    '''With a given list of possible configurations of cubes find the one that
    has the highest tota alignment score'''

    sums = []
    for config in configs:
        score_sum = 0
        for cube in config:
            score_sum += cube.align_score
        sums.append(score_sum)

    return configs[sums.index(max(sums))]



class Cube:

    def __init__(self, motif, align_score, pos1, pos2, read, quality_score = -1, second = None, morph = None):
        self.motif = motif
        self.align_score = align_score
        self.pos1 = pos1
        self.pos2 = pos2
        self.quality_score = quality_score
        self.read = read
        self.second = second
        self.morph = morph

    def __len__(self):
        return self.pos2 - self.pos1

    def __repr__(self):
        return f'{self.motif}, Score: {self.quality_score}, Range: {self.pos1} - {self.pos2}'

    def is_overlapping(self, c):
        # check if the boundaries of this cube overlap with the boundaries of another cube
        if (self.pos2 >= c.pos1 and self.pos1 <= c.pos1) or (self.pos1 <= c.pos2 and self.pos2 >= c.pos2):
            return True
        return False

    def calculate_qs(self):
        # calculate the quality score as alignment score squared divided by second-highest score
        score = 0
        second = ''
        for motif in motifs.keys():
            if not motif == self.motif:
                new_score = max(max(semiglobal_matrix(motif, self.read)[-1:][0])/len(motif),score)
                if new_score > score:
                    score = new_score
                    second = motif


        score = self.align_score ** 2 / score
        # score = self.score / score
        self.quality_score = score
        self.second = second

    def create_morph(self):

        if not self.second:
            self.calculate_qs()


        if self.quality_score < 1.2:
            # calculate a morphed motif as combination of the most and second-most-likely motif
            morph = ''
            for a,b in zip(translate(self.motif), translate(self.second)):
                if a == b:
                    morph += a
                else:
                    morph += '-'
    
            self.morph = morph
        else:
            self.morph = translate(self.motif)


def sort_into_closed_intervals(cubes):

    # figure out where unparsed regions are disjunct to avoid considering to many configurations later
    # sort the cubes and then find out where their intervals are disjunct

    # Sort intervals by the starting point
    cubes.sort(key = lambda x: x.pos1)

    # Iterate through the sorted intervals and check for gaps
    current_start, current_end = (cubes[0].pos1, cubes[0].pos2)
    current_id = 0

    cubes_sliced = []
    for i, cube in enumerate(cubes):
        start, end = (cube.pos1, cube.pos2)
        # If the start of the next interval is greater than the current end, there's a gap
        if start > current_end:
            cubes_sliced.append(cubes[current_id:i])
            current_id = i
        # Otherwise, merge the intervals by extending the current end
        current_end = max(current_end, end)

    cubes_sliced.append(cubes[current_id:])
    # If no gaps were found, the union is a single closed interval
    return cubes_sliced



# def parse_motifs(infile_path, outfile1_path, outfile2_path,  num_seqs):
def parse_motifs(infile_path, outfile_path, num_seqs):


    # dna sequences that mark the beginning and end of the HVD
    start = "GAAAGTGGTAAAAGACCCGGGAGC"
    end = "CATAAACACGGAGGTTATGACGAG"

    total = 0
    omitted = 0

    # ensure the output files are clear
    # with open(outfile1_path, 'w'), open(outfile2_path, 'w'):
    #     pass
    with open(outfile_path, 'w'):
        pass

    # with (open(infile_path, 'r') as infile, open(outfile1_path,'a') as outfile1, open(outfile2_path, 'a') as outfile2):
    with (open(infile_path, 'r') as infile, open(outfile_path, 'a') as outfile):

        outfile.write('read_id,parsed_read,quality_score,indel_score\n')
        # outfile.write('read_id,parsed_read,quality_score,perc_insertions,perc_deletions\n')
        n = 0

        for line in tqdm(infile):
            try:
                line = line.strip()
                if line.startswith('>'):  # New sequence header
                    # if the required number of sequences has been parsed, end the run
                    n += 1
                    if n > num_seqs:
                        break
                    read_id = line.strip()
                    continue

                # if n * 2 < 348: continue

                read = line.strip('\n')
                # find the start and end of the HVD
                start_pos, start_score = find_in_read(start, read)
                end_pos, end_score = find_in_read(end, read)

                # check if the HVD is found correctly
                has_good_scores = start_score  >= 0.6 and end_score >= 0.6

                if (end_pos - len(end) - start_pos) < 0 or not has_good_scores:
                    # if the start and end is not found properly, check the reverse compliment of the DNA
                    read = reverse_complement(read)

                    start_pos, start_score = find_in_read(start, read)
                    end_pos, end_score = find_in_read(end, read)

                    # check if the HVD is found correctly
                    has_good_scores = start_score >= 0.6 and end_score >= 0.6

                    # If the start and end of the HVD is still not found correctly, ommit the sequence
                    if (end_pos - len(end) - start_pos) < 0 or not has_good_scores:

                        omitted += 1
                        continue

                # Remove the conserved domains
                read = read[start_pos:(end_pos-len(end))]

                # Use the read to do the alignment, use a second string of the same length to store the replaced pieces
                # In the original read replace the parts that have been used with spaces

                # Create the second string
                # parsed_read = ' ' * len(read)

                # Find exact matches for the motifs where the motifs don't overlap
                read, parsed_read, quality_string = replace_exact_matches(read, motifs)

                # Do the alignment for each motif
                # forwards and backwards
                value_matrix_forwards = np.zeros((len(motifs.keys()), len(read) + 1))
                value_matrix_backwards = np.zeros((len(motifs.keys()), len(read) + 1))
                motifs_ordered = []
                cubes = []

                for i, motif in enumerate(motifs.keys()):
                    # store the order of the motifs to retrieve ids later
                    motifs_ordered.append(motif)
                    # Get the alignment matrix
                    alignment_matrix_forwards = semiglobal_matrix(motif, read)
                    alignment_matrix_backwards = semiglobal_matrix(motif[::-1], read[::-1])

                    # get normalized scores for the probability of this motif at each position
                    value_matrix_forwards[i, :] = alignment_matrix_forwards[-1, :] / len(
                            motif)
                    value_matrix_backwards[i, :] = alignment_matrix_backwards[-1, :] / len(
                            motif)

                    # Idea: Define Cube objects as ranges where a motif fits, the height is
                    # the normalized score of how well the motif fits the range and the
                    # width is the length of the range. then try to fit a configuration
                    # of cubes that maximizes the sum of scores

                    # get the scores for how well the motifs fit
                    # invert the backwards scores
                    scores_fwd = value_matrix_forwards[i,:]
                    scores_bwd = value_matrix_backwards[i,:][::-1]


                    tresh = 0.65
                    # delete the scores that are beneath a threshold 0.7
                    for j in range(len(scores_fwd)):
                        if scores_fwd[j] < tresh:
                            scores_fwd[j] = 0
                        if scores_bwd[j] < tresh:
                            scores_bwd[j] = 0


                    # get the first rev_back score and match with the first forwards score to build the cube
                    for j in range(len(scores_bwd)):
                        if scores_bwd[j] == 0:
                            continue

                        for k in range(j,len(scores_fwd)):
                            if scores_bwd[j] == scores_fwd[k] and (k - j) in range(10,25) :
                                cubes.append(Cube(motif, scores_bwd[j], j,k, read[j:k]))
                                break

                for j in range(len(cubes)):

                    # if two or more cubes are in the exact same position check if they have the same translation and if so, remove the one with the lower score
                    for k in range(j+1,len(cubes)):

                        if cubes[j].pos1 == cubes[k].pos1 and cubes[j].pos2 == cubes[k].pos2 and  \
                            translate(cubes[j].motif) == translate(cubes[k].motif) :

                            if cubes[k].align_score >= cubes[j].align_score:
                                cubes.pop(j)
                            else:
                                cubes.pop(k)

                        if k >= len(cubes)-1: break
                        if j >= len(cubes)-1: break

                    if j >= len(cubes) - 1: break

                    # if a cube is not overlapping with another cube, parse it.
                    # if all([not cubes[j].is_overlapping(c) or cubes[j] == c for c in cubes]):
                    #     m = (cubes[j].pos2 - cubes[j].pos1)
                    #     read = read[:cubes[j].pos1] + ' ' * m + read[cubes[j].pos2:]
                    #     parsed_read = parsed_read[:cubes[j].pos1] + ' ' * math.floor(m/3) \
                    #                   + translate(cubes[j].motif) + ' ' * (m - len(translate(cubes[j].motif))- math.floor(m/3)) + parsed_read[cubes[j].pos2:]
                    #

                # find out which cubes are actually overlapping and adress each intervall individually
                if cubes:

                    disjunct_interval_cubes = sort_into_closed_intervals(cubes)

                    # if one interval has to many possible configurations skip the reads (because its computationally expensive and likely to be a very messy read)
                    if any(len(interval) > 20 for interval in disjunct_interval_cubes):
                        omitted += 1
                        continue
                else:
                    disjunct_interval_cubes = []


                # else find the configuration that maximizes the average
                for interval in disjunct_interval_cubes:

                    # get all configurations and sort out those that are impossible because cubes are overlapping
                    cube_configs = get_possible_configurations(interval)

                    # take the configuration of motifs with the highest score
                    best_config = get_best_config(cube_configs)

                    # Calculate the quality scores and if the qs is low, create a morphed motif from the most likely and second most likely motif
                    for cube in best_config:
                        cube.calculate_qs()
                        cube.create_morph()

                    for cube in best_config:

                        # calculate the length of the cube
                        m = cube.pos2 - cube.pos1

                        # replace the parsed part of the read with nothing but keep the length of the read
                        read = read[:cube.pos1] + ' ' * m + read[cube.pos2:]
                        # insert the parsed section from the read into the parsed read
                        # parsed_read = parsed_read[:cube.pos1] + \
                        #               ' ' * math.floor(m/3) + \
                        #               translate(cube.motif) + \
                        #               ' ' * (m - len(translate(cube.motif)) - math.floor(m/3) ) + \
                        #               parsed_read[cube.pos2:]

                        parsed_read = parsed_read[:cube.pos1] + \
                                      ' ' * math.floor(m/3) + \
                                      cube.morph + \
                                      ' ' * (m - len(cube.morph) - math.floor(m/3) ) + \
                                      parsed_read[cube.pos2:]


                        # insert the quality score at the parsed position in the string for the quality scores
                        quality_string = quality_string[:cube.pos1] + \
                                        ' ' * math.floor(m/3) + \
                                        str(round(cube.quality_score,2)).ljust(len(translate(cube.motif))) + \
                                        ' ' * (m - len(translate(cube.motif)) - math.floor(m / 3)) + \
                                        quality_string[cube.pos2:]


                # Number of untranslated bases (insertions): # TODO make this work
                # num_insert = len(read.replace(' ', ''))
                # perc_insert = num_insert/len(read)
                # # Number of deletions
                # num_del = len(finalize_read(parsed_read).replace(' ', '')) * 3 - (len(read) - num_insert)
                # perc_del = num_del/len(read)




                outfile.write(read_id + ',')
                outfile.write(finalize_read(parsed_read)+ ',')
                outfile.write(finalize_read(quality_string)+ ',')
                outfile.write(f"{1-(len(read.replace(' ', ''))/len(read))}\n") # Ratio of nucleic acids that were parsed from the entire string
                # outfile.write(str(perc_insert) + ',') # TODO make this work
                # outfile.write(str(perc_del) + '\n')

            except Exception as inst:

                warnings.warn(f"{type(inst)}\n{inst.args}\n{inst}")
                omitted += 1

        print(f"{omitted}/{n}")

def main():

    # get the arguments from the command line
    parser = argparse.ArgumentParser(description="Count occurrences of specific strings in text files.")
    parser.add_argument('-i', '--input', required=True, help='File or folder with the sequences in fasta format')
    parser.add_argument('-o', '--aa_output_file', required=False, default='', help='Output file to write the parsed reads')
    parser.add_argument('-n', '--num_seqs', required=False, default=float('inf'), help='Number of sequences after which to automatically end the run' )

    args = parser.parse_args()

    ext = ''
    if not args.num_seqs == float('inf'):
        ext += '_n' + str(args.num_seqs)

    # If the input file is a folder run on every fasta file present
    if os.path.isdir(args.input):

        # Create folder for outputs
        output_folder = args.input + '_output'
        try:
            os.mkdir(output_folder)
        except:
            pass

        for file in os.listdir(args.input):
            filename = os.fsdecode(file)
            if filename.endswith(".fasta") or filename.endswith(".fa"):
                infile_path = args.input + '/' + filename
                # if not args.na_output_file == '':
                #     outfile_path = output_folder + '/' + os.path.splitext(filename)[0] +  '_' + args.na_output_file
                # else:
                outfile_path = output_folder + '/' + os.path.splitext(filename)[0] + '_' + 'parsed' + ext

                num_seqs = float(args.num_seqs)

                print(f"Parsing reads from file: {filename} ...")

                parse_motifs(infile_path, outfile_path, num_seqs)

                print(f"Saved parsed reads to {outfile_path}")


    # if the input is not a directory run for the given file
    else:

        infile_path = args.input
        print(f"Parsing reads from file: {infile_path} ...")

        # Create file names
        if not args.aa_output_file == '':
            outfile_path = args.aa_output_file
            # outfile2_path = args.aa_output_file
        else:
            outfile_path = os.path.splitext(infile_path)[0] + '_' + 'dna' + ext
            # outfile2_path = os.path.splitext(infile_path)[0] + '_' + 'aa' + ext

        print(f"Saving parsed reads to {outfile_path}")

        num_seqs = float(args.num_seqs)
        # parse_motifs(infile_path, outfile1_path, outfile2_path,
        #               num_seqs)
        parse_motifs(infile_path, outfile_path,
                     num_seqs)



if __name__== '__main__':
    main()

