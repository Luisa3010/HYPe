# Motif parser for HYP1 in G. Pallida

This is a script that uses the Smith-Waterman alignment algorithm to find motifs in the Hyper Variable Domain of the HYP1 gene. 
Run the script on a fasta file for reads of the HYP1 gene to get the deduced motifs as DNA and amino acid sequences.

## Quickstart: 

To get the deduced amino acid sequences run this on a fasta file containing DNA reads for HYP1.

    python3 grep_motifs_sw.py -i <<input.fasta>> 

## Run configurations

You can also specify the following flags:

-i <<input>>: can be a file or a directory. If it is a directory the script will run on all contained files with extensions .fa and .fasta
-o1 <<dna-output-file>> where to store the motif sequence as DNA, if not set, a default file name will be used
-o2 <<aa-output-file>> where to store the deduce motif sequence as amino acids, if not set, a default file name will be used
-n <<num-sequences>> the number of sequences from the file to be analyzed, e.g. -n 20 will analyze the first 20 sequences
-b <<allow-bad>> a boolean specifying whether sequences with ambiguous motifs should be included. Defaults to False. If set to True, ambiguous motifs are marked with 'X'. There will also be an additional 'X' appended in the last position to mark the entire sequence.


## Run on several files in parallel:

This script is designed to be compatible with GNU parallel processing. It can run on several files in parallel.
To use it, first run:

    sudo apt-get install parallel

The most straightforward way to run this in parallel is to place the grep_motifs_sw.py file in a directory containing all the fasta files you want to analyze.
Then open this location in the terminal and run:

    parallel -j <<number-of-jobs>> python3 grep_motifs_sw.py -i {1}  ::: *.fasta *.fa  


But you can also specify further configurations as before like this:

    parallel -j <<number-of-jobs>> python3 grep_motifs_sw.py -i {1} -n {2} -b {3} ::: <<input-file1>> <<input-file2>> ... ::: <<num-sequences1>> <<num-sequences2>> ... ::: True False 

All the different possible combinations of the inputs will be used. Note, that because of this setting the -o1 and -o2 flags will not work correctly in parallel.
Also, if you give directories as inputs for the parallel run, each directory will be given to one thread, so it's most of the time advised to use it on the individual files rather than directories.

### Example:

    parallel -j 10 python3 grep_motifs_sw.py -i {1} -n {2} -b {3} ::: heads*.fasta ::: 10 100 ::: True 

Will run 10 jobs in parallel on all files starting with 'heads' and ending in '.fasta'. Each of these will be run for the first 10 and the first 100 parameters with errors allowed.

    prallel -j 2 python3 grep_motifs_sw.py -i {1} -b {2} ::: file1.fasta file2.fasta ::: True False

This will run two jobs on the file1.fasta and file2.fasta once with errors allowed, and once with errors not allowed.






