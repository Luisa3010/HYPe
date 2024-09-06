###### calculate secondary strucures and energy for a given motif combination -
# made by Lukas Burhardt

import argparse
from tqdm import tqdm
import RNA
import csv

motifs = {
        "YERGGG1": "TATGAGCGCGGAGGCGGA",
        "YERGGG2": "TATGAGCGCGGAGGCGGG",
        "YERGGG3": "TATGAACGCGGAGGCGGA",
        "YERGGG4": "TATGAGCGTGGAGGCGGA",
        "YERGGG5": "TATGAACGCGGAGGCGGG",
        "YQRGGG1": "TATCAGCGCGGAGGCGGA",
        "YQRGGG2": "TATCAGCGCGGAGGCGGG",
        "SNRGGG1": "AGTAACCGCGGAGGCGGA",
        "SNRGGG2": "AGTAACCGCGGGGGCGGA",
        "SNRGGG3": "AGTAACCGCGGAGGCGGG",
        "SNRGGG4": "AGTAACCGCGGGGGCGGG",
        "SNRGGG5": "AGTAACCGTGGAGGCGGA",
        "SDRGD1": "AGTGACCGCGGAGAC",
        "SDRGD2": "AGTGACCGCGGAGAT",
        "SSRGD1": "AGTAGCCGCGGAGAC",
        "RDRGD1": "CGTGACCGCGGAGAC",
        "SDRGGG1": "AGTGACCGCGGAGGCGGA",
        "SDRGGG2": "AGCGACCGCGGAGGCGGA",
        "SDRGGG3": "AGTGACCGCGGAGGCGGG",
        "RDNKRG1": "CGTGACAATAAGCGCGGA",
        "REGGD1": "CGTGAAGGCGGAGAC",
        "RDRGGG1": "CGTGACCGCGGAGGCGGA",
        "SDRGE1": "AGTGACCGCGGAGAG",
        "GNRGGG1": "GGTAACCGCGGAGGCGGG",
        "GNRGGG2": "GGTAACCGCGGAGGCGGA",
        "GNRGGG3": "GGTAACCGCGGGGGCGGA",
        "RDDQRG1": "CGTGACGATCAGCGCGGA"
        # "AGATAACCGCGGAGGCGGA": "XSNRGGG"

}


def calculate_gibbs_energy(sequence):
    # Treat the sequence as DNA directly
    fc = RNA.fold_compound(sequence)
    # Get the MFE (minimum free energy) structure and its energy
    structure, mfe = fc.mfe()
    
    # Get the partition function and the free energy of the ensemble
    fc.pf()
    gfe = fc.eval_structure(structure)
    
    return structure, gfe

def calculate_energy_from_fasta(filename):
    # sequences = []
    with open(filename, "r") as f, open("heads_0_quad_energy_out.txt", "w") as out:
        for line in tqdm(f):
            quad = (line.strip(": 0\n").split(" "))[:4]
            for i in range(0,len(quad)):
                if quad[i] == "YERGGG6":
                    quad[i] = "YQRGGG1"
                if quad[i] == "YERGGG7":
                    quad[i] = "YQRGGG2"
            sequence_id = "".join(quad)
            sequence = "".join([motifs[motif] for motif in quad])
            structure, gfe = calculate_gibbs_energy(sequence)
            out.write(
                f"{sequence_id}\t{sequence}\t{0}\t{structure}\t{gfe}\n")
            # sequences.append({"id": sequence_id, "sequence": sequence, "count": 0})

def main():
    parser = argparse.ArgumentParser(
        description = "Count occurrences of specific strings in text files.")
    parser.add_argument('-i', '--input', required = True,
                        help = 'File or folder with the sequences in fasta format')


    args = parser.parse_args()

    # Input file containing multiple sequences
    input_file = args.input

    # Output files
    intermediate_output_file = "processed_output_with_sequences_and_count.txt"
    final_output_file = "sorted_output_with_sequences_and_count.txt"

    # Read sequences from the input file while maintaining order
    sequences = calculate_energy_from_fasta(input_file)

    # Open intermediate output file for writing


    print(f"Results saved to heads_0_quad_energy_out.txt")

if __name__ == "__main__":
    main()
#
# # Read and parse the intermediate output file
# data = []
# with open(intermediate_output_file, "r") as infile:
#     reader = csv.reader(infile, delimiter="\t")
#     headers = next(reader)
#     for row in reader:
#         if row:
#             # Parse Gibbs Free Energy (GFE) and store as a float for sorting
#             gfe = float(row[4].replace(" kcal/mol", ""))
#             data.append([row[0], row[1], row[2], row[3], gfe])
#
# # Sort data by Gibbs Free Energy (GFE)
# data_sorted = sorted(data, key=lambda x: x[4])
#
# # Write the sorted data to the final output file
# with open(final_output_file, "w") as outfile:
#     writer = csv.writer(outfile, delimiter="\t")
#     # Write headers
#     writer.writerow(["Sequence ID", "Sequence", "Count", "Secondary Structure", "Gibbs Free Energy (GFE)"])
#     # Write sorted data
#     for entry in data_sorted:
#         writer.writerow([entry[0], entry[1], entry[2], entry[3], f"{entry[4]:.2f} kcal/mol"])
#
# print(f"Sorted results saved to {final_output_file}")
