###### Create tables for each of the given files of parsed motif sequences
# how often each motif is followed by each other motif

import numpy as np
import pandas as pd
def read_counts(file_path):
    counts_by_files = {}
    counts = {}
    current_file = ''
    with open(file_path, 'r') as f:
        for line in f:
            if not current_file == '' and line[-8:-1] == ".fastq:":
                counts_by_files[current_file] = counts
            if line[-8:-1] == ".fastq:":
                current_file = (line.strip().split(': ')[0]).split('.')[0]
                counts = {}
            else:
                string, count = line.strip().split(': ')
                counts[string] = int(count)
        counts_by_files[current_file] = counts

    return counts_by_files


def get_motifs(string, all_motifs):
    for motif in all_motifs:
        if string[:len(motif)] == motif:
            motif1 = motif
            break
    for motif in all_motifs:
        if string[-len(motif):] == motif:
            motif2 = motif
            break

    return motif1, motif2


# get the counts and the motifs
infile = "/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT/combs_counts.txt"
motifs_file = "/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT/HYP1_motif_variants.txt"

out_path = "/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT"

with open(motifs_file) as file:
    motifs = [line.strip() for line in file]

counts_by_file = read_counts(infile)


for header, dict in counts_by_file.items():
    A = np.zeros((len(motifs), len(motifs)))
    df = pd.DataFrame(A, index = motifs, columns = motifs)
    for key, value in dict.items():
        motif1, motif2 = get_motifs(key, motifs)
        df[motif2].loc[motif1] = value

    df.to_csv(out_path + f"/{header}_counts.csv")



