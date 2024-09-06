###### Get stats how well different segmentations of the amino acid sequences into
# motifs work using segmentations from a file that is created using a kmer based
# approach


import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm


def load_sequences(filename):
    # This function loads sequences from a given filename
    with open(filename, 'r') as file:
        sequences = file.readlines()
    return [s.strip() for s in sequences]


def count_motifs(sequences, known_segments):
    # This function counts all segments using lists
    segment_counts = {}

    for seq in sequences:
        for segment in sorted(known_segments, key = len, reverse = True):
            while segment in seq:
                seq = seq.replace(segment, ' ', 1)
                if segment in segment_counts:
                    segment_counts[segment] += 1
                else:
                    segment_counts[segment] = 1
        known_segments.extend(seq.strip().split())
        for segment in sorted(known_segments, key = len, reverse = True):
            while segment in seq:
                seq = seq.replace(segment, ' ', 1)
                if segment in segment_counts:
                    segment_counts[segment] += 1
                else:
                    segment_counts[segment] = 1


    return segment_counts


def categorize_motifs(segment_counts, threshold = 2):
    # This function categorizes motifs based on their occurrence
    # First, sort the segment counts dictionary by values in descending order
    sorted_segments = sorted(segment_counts.items(), key = lambda x: x[1],
                             reverse = True)

    # Now split them into true and rogue motifs based on the threshold
    true_motifs = {seg: count for seg, count in sorted_segments if
                   count > threshold}
    rogue_motifs = {seg: count for seg, count in sorted_segments if
                    count <= threshold}

    return true_motifs, rogue_motifs


def create_plot(true_motifs, rogue_motifs, filename):
    # This function creates a plot for true and rogue motif counts
    fig, ax = plt.subplots()
    ax.bar(true_motifs.keys(), true_motifs.values(), color = 'b',
           label = 'True Motifs')
    ax.bar(rogue_motifs.keys(), rogue_motifs.values(), color = 'r',
           label = 'Rogue Motifs')
    plt.xticks(rotation=90)
    plt.tight_layout(pad= 2.0)
    ax.set_xlabel('Segments')
    ax.set_ylabel('Frequency')
    ax.set_title('Motif Frequency Analysis')
    ax.legend()
    plt.savefig(filename)
    plt.close()


def create_csv(data, filename):
    # This function creates a CSV file from the provided data
    df = pd.DataFrame(data)
    df.to_csv(filename, index = False)


# segmentations = [['KKYGS', 'KKYDS', 'EEKEEE', 'DDKEKE', 'DEKEEE', 'KKYGK', 'KKYGN', 'DEKEKE'],
# ['KKYGS', 'EEKEEE', 'DKEKE', 'DEKEEE', 'KKYGK', 'KKYDSD', 'KKYGN', 'DEKEKE'],
# ['DEKEE', 'EKKYGK', 'EEKEE', 'EKKYD', 'DEKEK', 'EKKYGN', 'KKYGN', 'EKKYGS', 'SDDKEK'],
# ['DEKEE', 'EKKYGK', 'EEKEE', 'DEKEK', 'DDKEK', 'EKKYGN', 'EKKYDS', 'KKYGN', 'EKKYGS'],
# ['DEKEE', 'EKKYGK', 'KKYGS', 'KKYDS', 'EEKEEE', 'EKKYGN', 'DDKEKE', 'KKYGN', 'DEKEKE'],
# ['DEKEE', 'EKKYGK', 'KKYGS', 'EEKEEE', 'DKEKE', 'EKKYGN', 'KKYDSD', 'KKYGN', 'DEKEKE'],
# ['DEKEE', 'EKKYGK', 'KKYGS', 'KKYDS', 'EEKEEE', 'DDKEKE', 'DEKEEE', 'KKYGN', 'DEKEKE'],
# ['DEKEE', 'EKKYGK', 'KKYGS', 'EEKEEE', 'DKEKE', 'DEKEEE', 'KKYDSD', 'KKYGN', 'DEKEKE'],
# ['EEKEE', 'EKKYD', 'DEKEK', 'EKKYGN', 'DEKEEE', 'KKYGK', 'KKYGN', 'EKKYGS', 'SDDKEK'],
# ['DEKEKE', 'EEKEE', 'EKKYD', 'DEKEK', 'DEKEEE', 'KKYGK', 'KKYGN', 'EKKYGS', 'SDDKEK'],
# ['EEKEE', 'DEKEK', 'EKKYDS', 'DDKEK', 'EKKYGN', 'DEKEEE', 'KKYGK', 'KKYGN', 'EKKYGS'],
# ['EEKEE', 'DEKEK', 'EKKYDS', 'DDKEK', 'DEKEEE', 'KKYGK', 'KKYGN', 'EKKYGS', 'DEKEKE'],
# ['EEKEE', 'KKYDS', 'DEKEK', 'DDKEK', 'DEKEEE', 'KKYGK', 'KKYGN', 'EKKYGS', 'DEKEKE'],
# ['KKYGS', 'EEKEE', 'KKYDS', 'DDKEK', 'DEKEEE', 'KKYGK', 'KKYGN', 'EKKYGS', 'DEKEKE'],
# ['KKYGS', 'EEKEE', 'KKYDS', 'DDKEKE', 'DEKEEE', 'KKYGK', 'KKYGN', 'EKKYGS', 'DEKEKE'],
# ['KKYGS', 'EEKEE', 'DKEKE', 'DEKEEE', 'KKYGK', 'KKYDSD', 'KKYGN', 'EKKYGS', 'DEKEKE'],
# ['KKYGS', 'EKKYD', 'EEKEEE', 'DEKEK', 'DEKEEE', 'KKYGK', 'SDDKEKE', 'KKYGN', 'DEKEKE'],
# ['KKYGS', 'EEKEEE', 'DEKEK', 'EKKYDS', 'DDKEKE', 'DEKEEE', 'KKYGK', 'KKYGN', 'DEKEKE'],
# ['KKYGS', 'EEKEEE', 'DEKEK', 'DKEKE', 'DEKEEE', 'KKYGK', 'KKYGN', 'EKKYDSD', 'DEKEKE'],
# ['KKYGS', 'KKYDS', 'EEKEEE', 'DDKEK', 'DEKEEE', 'KKYGK', 'KKYGN', 'EKKYGS', 'DEKEKE'],
# ['DEKEKE', 'KKYGS', 'EEKEEE', 'KYDSD', 'DKEKE', 'DEKEEE', 'KKYGK', 'KKYGN', 'DEKEKEK'],
# ['KKYGS', 'EKEKKYG', 'KYDSDDK', 'EEKEEE', 'DEKEEE', 'KKYGK', 'SDEKEKE', 'KKYGN', 'DEKEKEK'],
# ['KKYGS', 'KKYDS', 'DDKEKE', 'EKEEE', 'DEKEEE', 'KKYGK', 'KKYGN', 'KKYGKE', 'DEKEKE'],
# ['KKYGS', 'DKEKE', 'EKEEE', 'DEKEEE', 'KKYGK', 'KKYDSD', 'KKYGN', 'KKYGKE', 'DEKEKE'],
# ['KKYGS', 'KKYDS', 'DDKEKE', 'KKYGKD', 'EKEEE', 'DEKEEE', 'KKYGN', 'KKYGKE', 'DEKEKE'],
# ['KKYGS', 'DKEKE', 'KKYGKD', 'EKEEE', 'DEKEEE', 'KKYDSD', 'KKYGN', 'KKYGKE', 'DEKEKE'],
# ['EKEKE', 'KKYGS', 'KKYDS', 'KKYGSD', 'DDKEKE', 'KKYGKD', 'EKEEE', 'KKYGKE', 'KKYGND'],
# ['EKEKE', 'KKYGS', 'KKYGSD', 'DKEKE', 'KKYGKD', 'EKEEE', 'KKYDSD', 'KKYGKE', 'KKYGND'],
# ['DEKEE', 'KEKKY', 'DSDDKE', 'KEEKEE', 'EKKYGS', 'EEKKYGK', 'GSDEKEK', 'KKYGN', 'NDEKE', 'EKKYG'],
# ['DEKEE', 'KEKKY', 'KEEKEE', 'SDEKEK', 'EKKYGS', 'DSDDKEK', 'EEKKYGK', 'KKYGN', 'NDEKE', 'EKKYG'],
# ['DEKEE', 'KEKKY', 'KEEKEE', 'DEKEK', 'EKKYGS', 'DSDDKEK', 'EEKKYGK', 'KKYGN', 'NDEKE', 'EKKYG'],
# ['DEKEE', 'KEEKEE', 'KEKKYD', 'SDEKEK', 'EKKYGS', 'EEKKYGK', 'KKYGN', 'NDEKE', 'EKKYG', 'SDDKEK'],
# ['DEKEE', 'KEEKEE', 'KEKKYD', 'DEKEK', 'EKKYGS', 'EEKKYGK', 'KKYGN', 'NDEKE', 'EKKYG', 'SDDKEK'],
# ['DEKEE', 'KEEKEE', 'KEKKYDS', 'DDKEK', 'SDEKEK', 'EKKYGS', 'EEKKYGK', 'KKYGN', 'NDEKE', 'EKKYG'],
# ['DEKEE', 'KEEKEE', 'KEKKYDS', 'DDKEK', 'DEKEK', 'EKKYGS', 'EEKKYGK', 'KKYGN', 'NDEKE', 'EKKYG'],
# ['DEKEE', 'EKKYD', 'KDEKEE', 'KEEKEE', 'DEKEK', 'EKKYGN', 'KKYGN', 'EKKYGS', 'EKKYG', 'SDDKEK'],
# ['DEKEE', 'KDEKEE', 'KEEKEE', 'DEKEK', 'DDKEK', 'EKKYGN', 'EKKYDS', 'KKYGN', 'EKKYGS', 'EKKYG'],
# ['DEKEE', 'KKYGS', 'KKYDS', 'KDEKEE', 'EKKYGN', 'DDKEKE', 'KEEKEEE', 'KKYGN', 'EKKYG', 'DEKEKE'],
# ['DEKEE', 'KKYGS', 'KDEKEE', 'DKEKE', 'EKKYGN', 'KEEKEEE', 'KKYDSD', 'KKYGN', 'EKKYG', 'DEKEKE'],
# ['DEKEE', 'EKKYGK', 'EKKYD', 'KEEKEE', 'DEKEK', 'EKKYGN', 'KKYGN', 'EKKYGS', 'EKKYG', 'SDDKEK'],
# ['DEKEE', 'EKKYGK', 'KEEKEE', 'DEKEK', 'DDKEK', 'EKKYGN', 'EKKYDS', 'KKYGN', 'EKKYGS', 'EKKYG'],
# ['DEKEE', 'EKKYGK', 'KKYGS', 'KKYDS', 'EKKYGN', 'DDKEKE', 'KEEKEEE', 'KKYGN', 'EKKYG', 'DEKEKE'],
# ['DEKEE', 'EKKYGK', 'KKYGS', 'DKEKE', 'EKKYGN', 'KEEKEEE', 'KKYDSD', 'KKYGN', 'EKKYG', 'DEKEKE'],
# ['DEKEE', 'EKKYGK', 'EEKEE', 'EKKYD', 'KEKKYGS', 'DEKEK', 'SDDKE', 'EKKYGN', 'KKYGN', 'EKKYGS'],
# ['DEKEKE', 'DEKEE', 'EKKYGK', 'EEKEE', 'EKKYD', 'DEKEK', 'EKKYGN', 'KKYGN', 'EKKYGS', 'SDDKEK'],
# ['DEKEE', 'EKKYGK', 'EEKEE', 'EKKYD', 'KKYGS', 'DEKEK', 'EKKYGN', 'SDDKEKE', 'KKYGN', 'EKKYGS'],
# ['DEKEE', 'EKKYGK', 'EEKEE', 'DEKEK', 'DDKEK', 'EKKYGN', 'EKKYDS', 'KKYGN', 'EKKYGS', 'DEKEKE'],
# ['DEKEE', 'EKKYGK', 'EEKEE', 'KKYGS', 'DEKEK', 'EKKYGN', 'DDKEKE', 'EKKYDS', 'KKYGN', 'EKKYGS'],
# ['DEKEE', 'EKKYGK', 'EEKEE', 'KKYGS', 'EKKYDSD', 'DEKEK', 'DKEKE', 'EKKYGN', 'KKYGN', 'EKKYGS'],
# ['DEKEE', 'EKKYGK', 'EEKEE', 'KKYDS', 'DEKEK', 'DDKEK', 'EKKYGN', 'KKYGN', 'EKKYGS', 'DEKEKE'],
# ['DEKEE', 'EKKYGK', 'EEKEE', 'KKYDS', 'KKYGS', 'DDKEK', 'EKKYGN', 'KKYGN', 'EKKYGS', 'DEKEKE'],
# ['DEKEE', 'EKKYGK', 'EEKEE', 'KKYDS', 'KKYGS', 'EKKYGN', 'DDKEKE', 'KKYGN', 'EKKYGS', 'DEKEKE'],
# ['DEKEE', 'EKKYGK', 'EEKEE', 'KKYGS', 'DKEKE', 'EKKYGN', 'KKYDSD', 'KKYGN', 'EKKYGS', 'DEKEKE'],
# ['DEKEE', 'EKKYGK', 'EEKEE', 'KKYDSDD', 'KEKEKKY', 'EKKYGN', 'GSDEKEK', 'KKYGN', 'EKKYGS', 'DEKEKE'],
# ['DEKEE', 'EKKYGK', 'EEKEE', 'KEKEKKY', 'EKKYGN', 'KYDSDD', 'GSDEKEK', 'KKYGN', 'EKKYGS', 'DEKEKEK'],
# ['DEKEE', 'EKKYGK', 'EEKEE', 'KYDSDDK', 'EKKYGN', 'EKEKKY', 'GSDEKEK', 'KKYGN', 'EKKYGS', 'DEKEKEK'],
# ['DEKEE', 'EKKYGK', 'EEKEE', 'EKEKKYG', 'KYDSDDK', 'SDEKEK', 'EKKYGN', 'KKYGN', 'EKKYGS', 'DEKEKEK'],
# ['DEKEE', 'EKKYGK', 'EEKEE', 'EKKYD', 'DEKEK', 'EKKYGN', 'DEKEEE', 'KKYGN', 'EKKYGS', 'SDDKEK'],
# ['DEKEE', 'EKKYGK', 'EEKEE', 'DEKEK', 'EKKYDS', 'DDKEK', 'EKKYGN', 'DEKEEE', 'KKYGN', 'EKKYGS'],
# ['DEKEE', 'EKKYGK', 'KKYGS', 'EKKYD', 'EEKEEE', 'DEKEK', 'EKKYGN', 'KKYGN', 'EKKYGS', 'SDDKEK'],
# ['DEKEE', 'EKKYGK', 'KKYGS', 'EKKYD', 'EEKEEE', 'DEKEK', 'EKKYGN', 'SDDKEKE', 'KKYGN', 'EKKYGS'],
# ['DEKEE', 'EKKYGK', 'KKYGS', 'EKKYD', 'EEKEEE', 'DEKEK', 'EKKYGN', 'SDDKEKE', 'KKYGN', 'DEKEKE'],
# ['DEKEE', 'EKKYGK', 'KKYGS', 'EEKEEE', 'DEKEK', 'DDKEK', 'EKKYGN', 'EKKYDS', 'KKYGN', 'EKKYGS'],
# ['DEKEE', 'EKKYGK', 'KKYGS', 'EEKEEE', 'DEKEK', 'EKKYGN', 'DDKEKE', 'EKKYDS', 'KKYGN', 'EKKYGS'],
# ['DEKEE', 'EKKYGK', 'KKYGS', 'EEKEEE', 'DEKEK', 'EKKYGN', 'DDKEKE', 'EKKYDS', 'KKYGN', 'DEKEKE'],
# ['DEKEE', 'EKKYGK', 'KKYGS', 'EEKEEE', 'DEKEK', 'DKEKE', 'EKKYGN', 'EKKYGS', 'KKYGN', 'EKKYDSD'],
# ['DEKEE', 'EKKYGK', 'KKYGS', 'EEKEEE', 'DEKEK', 'DKEKE', 'EKKYGN', 'KKYGN', 'EKKYDSD', 'DEKEKE'],
# ['DEKEE', 'EKKYGK', 'KKYGS', 'KKYDS', 'EEKEEE', 'DDKEK', 'EKKYGN', 'KKYGN', 'EKKYGS', 'DEKEKE'],
# ['DEKEE', 'EKKYGK', 'KKYGS', 'KKYDS', 'EEKEEE', 'DEKEK', 'EKKYGN', 'DDKEKE', 'KKYGN', 'DEKEKE'],
# ['DEKEE', 'EKKYGK', 'KKYGS', 'EEKEEE', 'DKEKE', 'DEKEK', 'EKKYGN', 'KKYDSD', 'KKYGN', 'DEKEKE'],
# ['DEKEKE', 'DEKEE', 'EKKYGK', 'KKYGS', 'EEKEEE', 'KYDSD', 'DKEKE', 'EKKYGN', 'KKYGN', 'DEKEKEK'],
# ['DEKEE', 'EKKYGK', 'KKYGS', 'EKEKKYG', 'KYDSDDK', 'EEKEEE', 'EKKYGN', 'SDEKEKE', 'KKYGN', 'DEKEKEK'],
# ['DEKEE', 'EKKYGK', 'KKYGS', 'KKYDS', 'EEKEEE', 'EKKYGN', 'DDKEKE', 'DEKEEE', 'KKYGN', 'DEKEKE'],
# ['DEKEE', 'EKKYGK', 'KKYGS', 'EEKEEE', 'DKEKE', 'EKKYGN', 'DEKEEE', 'KKYDSD', 'KKYGN', 'DEKEKE'],
# ['DEKEE', 'EKKYGK', 'KKYGS', 'KKYDS', 'EKKYGKE', 'EKKYGN', 'DDKEKE', 'EKEEE', 'KKYGN', 'DEKEKE'],
# ['DEKEE', 'EKKYGK', 'KKYGS', 'DKEKE', 'EKKYGKE', 'EKKYGN', 'EKEEE', 'KKYDSD', 'KKYGN', 'DEKEKE'],
# ['DEKEE', 'KKYGS', 'KKYDS', 'EKKYGKD', 'EKKYGN', 'DDKEKE', 'EKEEE', 'KKYGN', 'KKYGKE', 'DEKEKE'],
# ['DEKEE', 'KKYGS', 'EKKYGKD', 'DKEKE', 'EKKYGN', 'EKEEE', 'KKYDSD', 'KKYGN', 'KKYGKE', 'DEKEKE'],
# ['DEKEE', 'EEKEE', 'EKKYD', 'DEKEK', 'EKKYGN', 'DEKEEE', 'KKYGK', 'KKYGN', 'EKKYGS', 'SDDKEK'],
# ['DEKEE', 'EEKEE', 'DEKEK', 'EKKYDS', 'DDKEK', 'EKKYGN', 'DEKEEE', 'KKYGK', 'KKYGN', 'EKKYGS'],
# ['DEKEE', 'KKYGS', 'KKYDS', 'EEKEEE', 'EKKYGN', 'DDKEKE', 'DEKEEE', 'KKYGK', 'KKYGN', 'DEKEKE'],
# ['DEKEE', 'KKYGS', 'EEKEEE', 'DKEKE', 'EKKYGN', 'DEKEEE', 'KKYGK', 'KKYDSD', 'KKYGN', 'DEKEKE'],
# ['DEKEE', 'KKYGS', 'KKYDS', 'KDEKEE', 'DDKEKE', 'DEKEEE', 'KEEKEEE', 'KKYGN', 'EKKYG', 'DEKEKE'],
# ['DEKEE', 'KKYGS', 'KDEKEE', 'DKEKE', 'DEKEEE', 'KEEKEEE', 'KKYDSD', 'KKYGN', 'EKKYG', 'DEKEKE'],
# ['DEKEE', 'EKKYGK', 'KKYGS', 'KKYDS', 'DDKEKE', 'DEKEEE', 'KEEKEEE', 'KKYGN', 'EKKYG', 'DEKEKE'],
# ['DEKEE', 'EKKYGK', 'KKYGS', 'DKEKE', 'DEKEEE', 'KEEKEEE', 'KKYDSD', 'KKYGN', 'EKKYG', 'DEKEKE'],
# ['DEKEKE', 'DEKEE', 'EKKYGK', 'EEKEE', 'EKKYD', 'DEKEK', 'DEKEEE', 'KKYGN', 'EKKYGS', 'SDDKEK'],
# ['DEKEE', 'EKKYGK', 'EEKEE', 'DEKEK', 'EKKYDS', 'DDKEK', 'DEKEEE', 'KKYGN', 'EKKYGS', 'DEKEKE'],
# ['DEKEE', 'EKKYGK', 'EEKEE', 'KKYDS', 'DEKEK', 'DDKEK', 'DEKEEE', 'KKYGN', 'EKKYGS', 'DEKEKE'],
# ['DEKEE', 'EKKYGK', 'EEKEE', 'KKYDS', 'KKYGS', 'DDKEK', 'DEKEEE', 'KKYGN', 'EKKYGS', 'DEKEKE'],
# ['DEKEE', 'EKKYGK', 'EEKEE', 'KKYDS', 'KKYGS', 'DDKEKE', 'DEKEEE', 'KKYGN', 'EKKYGS', 'DEKEKE'],
# ['DEKEE', 'EKKYGK', 'EEKEE', 'KKYGS', 'DKEKE', 'DEKEEE', 'KKYDSD', 'KKYGN', 'EKKYGS', 'DEKEKE'],
# ['DEKEE', 'EKKYGK', 'KKYGS', 'EKKYD', 'EEKEEE', 'DEKEK', 'DEKEEE', 'SDDKEKE', 'KKYGN', 'DEKEKE'],
# ['DEKEE', 'EKKYGK', 'KKYGS', 'EEKEEE', 'DEKEK', 'EKKYDS', 'DDKEKE', 'DEKEEE', 'KKYGN', 'DEKEKE'],
# ['DEKEE', 'EKKYGK', 'KKYGS', 'EEKEEE', 'DEKEK', 'DKEKE', 'DEKEEE', 'KKYGN', 'EKKYDSD', 'DEKEKE'],
# ['DEKEE', 'EKKYGK', 'KKYGS', 'KKYDS', 'EEKEEE', 'DDKEK', 'DEKEEE', 'KKYGN', 'EKKYGS', 'DEKEKE'],
# ['DEKEKE', 'DEKEE', 'EKKYGK', 'KKYGS', 'EEKEEE', 'KYDSD', 'DKEKE', 'DEKEEE', 'KKYGN', 'DEKEKEK'],
# ['DEKEE', 'EKKYGK', 'KKYGS', 'EKEKKYG', 'KYDSDDK', 'EEKEEE', 'DEKEEE', 'SDEKEKE', 'KKYGN', 'DEKEKEK'],
# ['DEKEE', 'EKKYGK', 'KKYGS', 'KKYDS', 'EKKYGKE', 'DDKEKE', 'EKEEE', 'DEKEEE', 'KKYGN', 'DEKEKE'],
# ['DEKEE', 'EKKYGK', 'KKYGS', 'DKEKE', 'EKKYGKE', 'EKEEE', 'DEKEEE', 'KKYDSD', 'KKYGN', 'DEKEKE']]

segmentations = []

dir = '/home/luisa/Documents/Work/Uni Cambridge/Data/HYP3_sanger/'

with open(dir + 'kmers', 'r') as file:
    for line in file:
        segmentations.append(line.strip().split(','))

sequences = load_sequences(dir + 'HYP3_protein_repeat_only_cleaned_stripped_fasta.txt')

motif_combinations_found = list(range(len(segmentations)))
true_motifs = list(range(len(segmentations)))
rogue_motifs = list(range(len(segmentations)))


for i in tqdm(range(0, len(segmentations))):

    motif_combinations_found[i] = count_motifs(sequences, segmentations[i])
    # Categorize segments
    true_motifs[i], rogue_motifs[i] = categorize_motifs(motif_combinations_found[i])

# Analyze results
has_short_motifs = [False for i in range(0,len(segmentations))]
for i in range(0, len(segmentations)):
    for motif in true_motifs[i]:
        if len(motif) == 1 or len(motif) == 2 or len(motif) == 3:
            has_short_motifs[i] = True

data = {
        'Variant ID': range(1, len(segmentations) + 1),
        'num_kmers': [len(motif_combinations_found[i]) for i in range(0, len(segmentations))],
        'num_true_motifs': [len(true_motifs[i]) for i in range(0, len(segmentations))],
        'num_rogue_motifs': [len(rogue_motifs[i]) for i in range(0, len(segmentations))],
        'has_short_motifs': has_short_motifs,
        'kmers': [[true_motifs[i],rogue_motifs[i]] for i in range(0, len(segmentations))]

}


create_csv(data, dir + 'segmentation_analysis.csv')
for i in range(0, len(segmentations)):
    if len(true_motifs[i]) <= 16 and not has_short_motifs[i] :
        create_plot(true_motifs[i], rogue_motifs[i], dir + f'plots/segmentation{i+1}')



