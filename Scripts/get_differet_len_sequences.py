##### get a feeling what proportion of the reads in the plasmid mix with 0
# dilution are actually what they should be and look at the reads that are different

ref_seq = "YERGGG SNRGGG RDRGD YERGGG SNRGGG SDRGE YERGGG SNRGGG SDRGD YERGGG SNRGGG RDRGD YERGGG"

fasta_file = "/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/semiglobal_test_run/plasmid_mix1_rep1_dil0in1000.phred15.400to1400bp.alignscore2000.fulllength_aa_n2843"
with open(fasta_file, 'r') as file:
    for line in file:
        if line.startswith('>'):
            id = line.strip("\n")
        else:
            line = line.strip('\n')
            if not len(line) == len(ref_seq):
                print(id)
                print(line)


