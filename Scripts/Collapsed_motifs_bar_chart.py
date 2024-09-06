##### Script to explore if reducing the motifs to a few clusters helps with getting
# a higher amount of the correct sequences in the plasmid 0 control

import os
from collections import Counter
from Bio import SeqIO
import matplotlib.pyplot as plt
import re

def natural_sort_key(s):
    return [float(text) if re.match(r'^\d*\.?\d+$', text) else text.lower() for text in re.split(r'(\d+\.\d+|\d+)', s)]


def read_fasta(file):
    sequences = []
    for record in SeqIO.parse(file, "fasta"):
        sequences.append(str(record.seq))
    return sequences

def count_sequences(sequences):
    total_count = len(sequences)
    sequence_counts = Counter(sequences)
    sequence_proportions = {seq: count / total_count for seq, count in sequence_counts.items()}
    return sequence_proportions



def get_top_sequences(sequence_proportions, top_n=5):
    sorted_sequences = sorted(sequence_proportions.items(), key=lambda item: item[1], reverse=True)
    return sorted_sequences[:top_n]


def plot_bar_chart(file, top_sequences):
    sequences, counts = zip(*top_sequences)
    truncated_sequences = [f"{i + 1}: {seq[:25]}..." if len(seq) > 25 else f"{i+1}: {seq[:25]}" for i, seq in enumerate(sequences)]
    plt.figure(figsize=(10, 6))
    plt.bar(truncated_sequences, counts, color='blue')
    plt.xlabel('Sequences')
    plt.ylabel('Proportion')
    plt.title(f'Top {len(sequences)} sequences in {os.path.basename(file)}')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()



def process_fasta_files(files):
    
    fig, ax = plt.subplots(2,len(files), constrained_layout = True)
    fig.suptitle("Top 5 Sequences in Mix1Rep1 original vs collapsed motifs")

    
    for i, file in enumerate(files):
        print(f'Processing file: {file}')
        sequences = read_fasta(file)
        sequence_counts = count_sequences(sequences)
        top_sequences = get_top_sequences(sequence_counts)
        
        
        # plot_bar_chart(file, top_sequences)
        seqs, counts = zip(*top_sequences)
        truncated_sequences = [f"{i + 1}: {seq[:15]}..." if len(
            seq) > 15 else f"{i + 1}: {seq[:15]}" for i, seq in
                               enumerate(seqs)]
        ax[0,i].bar(truncated_sequences, counts, color = 'blue')
        # ax[0,i].set_xlabel('Sequences')
        ax[0,i].set_ylabel('Proportion')
        ax[0,i].set_title(f'{os.path.basename(file)[21:-(80-24)]}in1000')
        ax[0,i].set_xlim([0,1])
        ax[0,i].set_xticks(ticks = [0,1,2,3,4],labels = truncated_sequences, rotation = 90)

        
        for seq, count in top_sequences:
            print(f'Sequence: {seq}, Count: {count}')
        print("\n")

        collapsed_sequences = []
        for seq in sequences:
            for motif in collapsed_motifs.keys():
                seq = seq.replace(motif, collapsed_motifs[motif])
            collapsed_sequences.append(seq)


        collapsed_sequence_counts = count_sequences(collapsed_sequences)
        top_collapsed_sequences = get_top_sequences(collapsed_sequence_counts)
        # plot_bar_chart(file, top_collapsed_sequences)

        collapsed_seqs, counts = zip(*top_collapsed_sequences)
        truncated_collapsed_sequences = [f"{i + 1}: {seq[:15]}..." if len(
                seq) > 15 else f"{i + 1}: {seq[:15]}" for i, seq in
                               enumerate(collapsed_seqs)]
        ax[1,i].bar(truncated_collapsed_sequences, counts, color = 'blue')
        # ax[1,i].set_xlabel('Sequences')
        ax[1,i].set_ylabel('Proportion')
        # ax[1,i].set_title(
        #    f'{os.path.basename(file)[21:-(80-24)]}in1000')
        ax[1, i].set_xlim([0, 1])
        ax[1,i].set_xticks(ticks = [0, 1, 2, 3, 4],
                            labels = truncated_collapsed_sequences, rotation = 90)

        for seq, count in top_collapsed_sequences:
            print(f'Sequence: {seq}, Count: {count}')
        print("\n")

    plt.xticks(rotation = 90)
    # plt.tight_layout()
    plt.show()


collapsed_motifs = {
        'YERGGG': 'YXRGGG',
        'YQRGGG': 'YXRGGG',
        'SNRGGG': 'SXRGGG',
        'SDRGGG': 'SXRGGG',
        'GNRGGG': 'SXRGGG',
        'SDRGD': 'XDRGD',
        'RDRGD': 'XDRGD',
        'SDRGE': 'XDRGD',
        'SSRGD': 'XDRGD',
        'RDNKRG': 'RDXXRG',
        'RDDQRG': 'RDXXRG',
}


# Load files
dir = "/home/luisa/Documents/Work/Uni_Cambridge/Data/smith_waterman_motifs/plasmid/half1_output"
fasta_files = []
for file in os.listdir(dir):
    filename = os.fsdecode(file)
    if filename.endswith('AA'):
        fasta_files.append(dir + '/' + filename)

fasta_files = sorted(fasta_files, key = natural_sort_key)

# Process each FASTA file
process_fasta_files(fasta_files)
