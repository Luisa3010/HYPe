### A helper script for the counts triples of motifs in heads in tails to
# aggregate the counts from all the individual reps into one single file

# Define the input and output files
path = "/home/luisa/Documents/Work/Uni_Cambridge/grep_motifs"
input_files = [f'{path}/heads{i}_3motifs_out.txt' for i in range(1, 11)]
output_file = path + '/heads_3motifs_out.txt'

# Initialize a dictionary to store the aggregated counts
counts = {}

# Iterate over each input file
for file_name in input_files:
    with open(file_name, 'r') as file:
        for line in file:
            fragment, count = line.strip().split(': ')
            count = int(count)
            if fragment in counts:
                counts[fragment] += count
            else:
                counts[fragment] = count

# Write the aggregated counts to the output file
with open(output_file, 'w') as file:
    for fragment, count in sorted(counts.items()):
        file.write(f'{fragment}: {count}\n')

print(f'Aggregated counts have been written to {output_file}')
