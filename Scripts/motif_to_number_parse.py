###### for a file of motif HVDs for  HYP1 G. Pallida reads transform the motifs
# on DNA base to the IDs (amino acid sequence + number)


import argparse

def parse_line(line):
    motif_variations = {
            "TATGAGCGCGGAGGCGGA": "YERGGG1",
            "TATGAGCGCGGAGGCGGG": "YERGGG2",
            "TATGAACGCGGAGGCGGA": "YERGGG3",
            "TATGAGCGTGGAGGCGGA": "YERGGG4",
            "TATGAACGCGGAGGCGGG": "YERGGG5",
            "AGTAACCGCGGAGGCGGA": "SNRGGG1",
            "AGTAACCGCGGGGGCGGA": "SNRGGG2",
            "AGTAACCGCGGAGGCGGG": "SNRGGG3",
            "AGTAACCGCGGGGGCGGG": "SNRGGG4",
            "AGTAACCGTGGAGGCGGA": "SNRGGG5",
            "CGTGACCGCGGAGAC": "RDRGD1",
            "AGTGACCGCGGAGAC": "SDRGD1",
            "AGTGACCGCGGAGAT": "SCRGD2",
            "AGTGACCGCGGAGGCGGA": "SDRGGG1",
            "AGCGACCGCGGAGGCGGA": "SDRGGG2",
            "AGTGACCGCGGAGGCGGG": "SDRGGG3",
            "CGTGACAATAAGCGCGGA": "RDNKRG1",
            "CGTGAAGGCGGAGAC": "REGGD1",
            "CGTGACCGCGGAGGCGGA": "RDRGGG1",
            "AGTGACCGCGGAGAG": "SDRGE1",
            "GGTAACCGCGGAGGCGGG": "GNRGGG1",
            "GGTAACCGCGGAGGCGGA": "GNRGGG2",
            "GGTAACCGCGGGGGCGGA": "GNRGGG3"
    }

    motifs = line.split(" ")
    parsed_motifs = [motif_variations[motifs[i]] for i in range(len(motifs))]
    return ' '.join(parsed_motifs)

def main():

    parser = argparse.ArgumentParser(
        description = "Count occurrences of specific strings in text files.")
    parser.add_argument('-i', '--input_file', required = True,
                        help = 'File with the dna sequences in fasta format')
    parser.add_argument('-o', '--output_file', required = True,
                        help = 'Output file to write the names with numbers of the motifs')

    args = parser.parse_args()

    infile_path = args.input_file
    outfile_path = args.output_file

    # outfile_path = "tmp"
    # infile_path = "/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/na_TEST.fasta"

    with open(outfile_path, 'w') as outfile:
        pass

    with open(infile_path, 'r') as infile, open(outfile_path, 'a') as outfile:
        current_seq_id = None
        sequence = ''

        for line in infile:

            line = line.strip()
            if line.startswith('>'):  # New sequence header
                # Extract the protein ID or some identifier from the header
                current_seq_id = line  # Takes the line as ID
                outfile.write(line + "\n")
                sequence = ''  # Reset sequence
            else:  # Append sequence without newlines
                parsed_line = parse_line(line)
                outfile.write(parsed_line + "\n")






if __name__ == "__main__":
    main()