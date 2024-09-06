####### Seperarete an amino acid sequence in all possible combinations of kmers
# using brute force, you can specify what size these kmers a re supposed to have


import time

def all_ways_to_segment_sequence(sequence, min_length, max_length):

    def generate_partitions(sequence, current_partition):
        if len(sequence) == 0:
            all_partitions.append(current_partition[:])
            return

        if len(list(set(current_partition))) > 14:
            return

        for i in range(min_length, min(max_length + 1, len(sequence) + 1)):
            current_partition.append(sequence[:i])
            generate_partitions(sequence[i:], current_partition)
            current_partition.pop()

    all_partitions = []
    generate_partitions(sequence, [])
    return all_partitions

start = time.time()

#sequence = "KKYGNDEKEEEKKYGKDEKEEEKKYGNDGKEKEKKYGSDENEKEE"
#sequence = "KKYGNDEKEEEKKYGNDEKEEEKKYGKDEKEEEKKYGKEEKEEEKKYGSDEKEEEKKYGNDEKEKEKKYDSDDKEKEKKYGSDEKEKEKKYG"
#sequence = "KKYGNDEKEEEKKYGDDEKEEEKKYGKDEKEEEKKYGKEEKEEEKKYGSDEKEEEKKYGNDEKEKEKKYDSDDKEKEKKYGSDEKEKEKKYGNDEKEKEKKYDSDDKEKEKKYGSDEKEKEKKYGSDENEKEE"
sequence = "KKYGNDEKEEEKKYGKDEKEEEKKYGKDEKEEEKKYGNEKEEEKKYGKDEKEEEKKYGSDEKEKKYGNDEKEEEKKYGKDEKEEEKKYGKDEKEEEKKYGSDEKEKEKKYGSDDKEEEKKYGNDEKEKEKKYGSDQKEEEKKYGKDEKEEEKKYGNDGKEKEKKYGSDENEKEE"
min_length = 5
max_length = 7
all_segmentations = all_ways_to_segment_sequence(sequence, min_length, max_length)

# Remove duplicates (if any) and sort by number of segments
unique_segmentations = []
seen = set()
for seg in all_segmentations:
    seg_tuple = tuple(sorted(list(set(seg))))
    if seg_tuple not in seen:
        seen.add(seg_tuple)
        unique_segmentations.append(list(set(seg)))

sorted_segmentations = sorted(unique_segmentations, key=len)

# Print top 10 segmentations with the fewest segments
for segmentation in sorted_segmentations[:100]:
    print(segmentation)

current_time  = time.strftime("%Y-%m-%d - %H-%M", time.localtime())


with open(f'/home/luisa/Documents/Work/Uni Cambridge/Data/kmers{current_time}','w') as file:
    for segmentation in sorted_segmentations:
        for i in range(0, len(segmentation)):
            file.write(segmentation[i])
            if not i == len(segmentation) -1:
                file.write(',')
        file.write('\n')



end = time.time()
print("Runtime:", end - start, "s")



# take the best segmentation
# delete the known segments from the other sequences
# count how many times each of them appear
# go over the remaining sequence and store all the remaining segments in a list
# count how many times each of those appears
# count how many true motifs and rogue there are and store it in a file
# create a plot of the counts and store it
# create table:
#   Variant ID  kmers               num_kmers    num_true_motifs num_rogue_motifs
#   1           ['ABC', 'DEF' ...]  5            3               2


