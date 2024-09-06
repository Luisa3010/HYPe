###### finds existing four motif combinations for G. Pallida HYP1

import argparse
def main():

    parser = argparse.ArgumentParser(
        description = "Count occurrences of specific strings in text files.")
    parser.add_argument('-i', '--input', required = True,
                        help = 'File or folder with the sequences in fasta format')
    parser.add_argument('-o1', '--non_zeros_output_file', required = True,
                        help = 'Output file to write the quadruples with count >0')
    parser.add_argument('-o2', '--zeros_output_file', required = True,
                        help = 'Output file to write the quadruples with count = 0')

    motifs_dict = {
            "TATGAGCGCGGAGGCGGA": "YERGGG1",
            "TATGAGCGCGGAGGCGGG": "YERGGG2",
            "TATGAACGCGGAGGCGGA": "YERGGG3",
            "TATGAGCGTGGAGGCGGA": "YERGGG4",
            "TATGAACGCGGAGGCGGG": "YERGGG5",
            "TATCAGCGCGGAGGCGGA": "YQRGGG1",
            "TATCAGCGCGGAGGCGGG": "YQRGGG2",
            "AGTAACCGCGGAGGCGGA": "SNRGGG1",
            "AGTAACCGCGGGGGCGGA": "SNRGGG2",
            "AGTAACCGCGGAGGCGGG": "SNRGGG3",
            "AGTAACCGCGGGGGCGGG": "SNRGGG4",
            "AGTAACCGTGGAGGCGGA": "SNRGGG5",
            "AGTGACCGCGGAGAC": "SDRGD1",
            "AGTGACCGCGGAGAT": "SDRGD2",
            "AGTAGCCGCGGAGAC": "SSRGD1",
            "CGTGACCGCGGAGAC": "RDRGD1",
            "AGTGACCGCGGAGGCGGA": "SDRGGG1",
            "AGCGACCGCGGAGGCGGA": "SDRGGG2",
            "AGTGACCGCGGAGGCGGG": "SDRGGG3",
            "CGTGACAATAAGCGCGGA": "RDNKRG1",
            "CGTGAAGGCGGAGAC": "REGGD1",
            "CGTGACCGCGGAGGCGGA": "RDRGGG1",
            "AGTGACCGCGGAGAG": "SDRGE1",
            "GGTAACCGCGGAGGCGGG": "GNRGGG1",
            "GGTAACCGCGGAGGCGGA": "GNRGGG2",
            "GGTAACCGCGGGGGCGGA": "GNRGGG3",
            "CGTGACGATCAGCGCGGA": "RDDQRG1",
            # "AGATAACCGCGGAGGCGGA": "XSNRGGG"

    }


    args = parser.parse_args()


    motifs_file = args.input

    # create a dictionary with all of the possible quadruples
    quadruples = {}
    for m1 in motifs_dict.keys():
        for m2 in motifs_dict.keys():
            for m3 in motifs_dict.keys():
                for m4 in motifs_dict.keys():
                    quadruples[' '.join([m1,m2,m3,m4])] = 0


    j = 0
    with open(motifs_file, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                j += 1
                if (j % 2000 == 0):
                    print(f"{j} Sequences done")
                motifs = line.split(' ')

                for i in range(len(motifs) -4):
                    quad = ' '.join(motifs[i:i+4])
                    quadruples[quad] += 1


    non_zeros_out_file = args.non_zeros_output_file
    zeros_out_file = args.zeros_output_file

    with open(non_zeros_out_file, 'w') as file:
        pass

    with open(zeros_out_file, 'w') as file:
        pass

    with open(non_zeros_out_file, 'a') as nz_file, open(zeros_out_file, 'a') as z_file:
        #
        # for key, value in quadruples.items():
        #     if value == 0:
        #         z_file.write(f"{key}: {value}\n")
        #     else:
        #         nz_file.write(f"{key}: {value}\n")


        for key, value in quadruples.items():
            motifs_split = key.split(' ')
            translated_motifs = []
            for motif in motifs_split:
                translated_motifs.append(motifs_dict[motif])
            trans_key = ' '.join(translated_motifs)
            if value == 0:
                z_file.write(f"{trans_key}: {value}\n")
            else:
                nz_file.write(f"{trans_key}: {value}\n")


if __name__ == "__main__":
    main()