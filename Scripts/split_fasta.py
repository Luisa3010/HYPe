##### split one large fasta file into several smaller files


import argparse

def split_fasta(input_file, num_files):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    num_sequences = len(lines) // 2
    sequences_per_file = num_sequences // num_files
    remainder = num_sequences % num_files

    output_files = [open(f"{input_file}_{i+1}.fasta", 'w') for i in range(num_files)]

    seq_index = 0
    for i in range(num_files):
        for j in range(sequences_per_file):
            output_files[i].write(lines[seq_index * 2])
            output_files[i].write(lines[seq_index * 2 + 1])
            seq_index += 1

        if i < remainder:
            output_files[i].write(lines[seq_index * 2])
            output_files[i].write(lines[seq_index * 2 + 1])
            seq_index += 1

    for f in output_files:
        f.close()

def main():
    parser = argparse.ArgumentParser(description='Split a FASTA file into multiple output files.')
    parser.add_argument('-i', '--inputfile', type=str, required=True, help='Input FASTA file')
    parser.add_argument('-n', '--numfiles', type=int, required=True, help='Number of output files')

    args = parser.parse_args()

    split_fasta(args.inputfile, args.numfiles)

if __name__ == '__main__':
    main()
