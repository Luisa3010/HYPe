##### Explore how much the positions of the "YERGGGs" differ to gauge how big
# or small potential cut-outs might be

def read_fasta(file_path):
    with open(file_path, 'r') as file:
        sequences = {}
        current_seq_id = None
        sequence = ''

        for line in file:
            line = line.strip()
            if line.startswith('>'):  # New sequence header
                if current_seq_id:
                    sequences[current_seq_id] = sequence
                # Extract the protein ID or some identifier from the header
                current_seq_id = line  # Takes the line as ID
                sequence = ''  # Reset sequence
            else:
                sequence += line.strip()  # Append sequence without newlines

        # Add the last sequence to the dictionary
        if current_seq_id:
            sequences[current_seq_id] = sequence

    return sequences



letter = "Y"

sequences = read_fasta('/home/luisa/Documents/Work/Uni Cambridge/Data/HYP1_protein_repeat_only_fasta.txt')
positions = []

for s, sequence in sequences.items():
    for index, character in enumerate(sequence):
        if character == letter:
            positions.append(index + 1)

print(sorted(positions))