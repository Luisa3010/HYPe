##### use a fuzzy matching approach to get the HVD from a G. Pallida HYP3 read



def is_close_match(s1, s2, max_diff = 1):
    """Check if two strings are close matches, allowing for a specified number of differences."""
    if len(s1) != len(s2):
        return False
    differences = sum(1 for a, b in zip(s1, s2) if a != b)
    return differences <= max_diff


def find_close_matches(text, target):
    """Find the first and last occurrence of substrings that match the target with at most one difference."""
    target_length = len(target)
    first_occurrence = -1
    last_occurrence = -1

    for i in range(len(text) - target_length + 1):
        substring = text[i:i + target_length]
        if is_close_match(substring, target):
            if first_occurrence == -1:
                first_occurrence = i
            last_occurrence = i

    return first_occurrence, last_occurrence

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
    while not len(seq) % 3 == 0:
        seq = seq +"A"
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        if codon in table.keys():
            protein += table[codon]
        else:
            return protein
    return protein
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


def cut_sequence(seq):
    # translate sequence
    translated_seq = translate(seq)

    try:
        # find index of first KKYG
        start, a = find_close_matches(translated_seq,"KKYG")
        # find index of SDENEKEE
        a,end = find_close_matches(translated_seq,"SDENEKEE")
        # cut DNA there
        seq = seq[start * 3 :(end + 8) *3 ]
    except:
        print("Warning: Couldn't cut a sequence")
    return seq

sequences = read_fasta("/home/luisa/Documents/Work/Uni Cambridge/Data/HYP3_sanger/HYP3_nucleotide_fasta")
cut_sequences = {}

for seq_id, sequence in sequences.items():
    cut_sequences[seq_id] = cut_sequence(sequence)



with open("tmp_output.txt", 'w') as file:
    for seq_id, sequence in cut_sequences.items():
        file.write(seq_id)
        file.write("\n")
        file.write(sequence)
        file.write("\n")