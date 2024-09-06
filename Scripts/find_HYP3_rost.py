###### find the HVDs for the reads that rouhly align with the HYP3 locus in G. Rostochiensis


from semiglobal import semiglobal_matrix
from grep_motifs_semiglobal import find_in_seqs, translate, complement_dna
import numpy as np
from tqdm import tqdm

# load all of the reads
path = '/home/luisa/Documents/Work/Uni_Cambridge/Data/Rostochiensis/PACBIO/L19/SRR1356039x_inbred_L19_reps123_mapq30_primaryalignonly_hypsonly_merged.fasta'
reads = []
read_ids = []
with open(path, 'r') as file:
    for line in file:
        if not line.startswith('>'):
            reads.append(line.strip('\n'))
        else:
            read_ids.append(line.strip('\n'))




# Define a consensus ssequence for the for the HVD
consensus_rev = 'catttcctGTCCATTTATTTAAACACAAACCTGCGCACAGCACGGAGACGATGCTCATACTTTTCTCCCTTGCTTTCCTTTTCATCGCTACCGTACTTCTCCTTTTCCTTTTCATCGTTACCGTACTTCTCCTTTTCCTTTTCATCGCTACCGTACTTCTCCTTTTCCTTTTCATCGTTACCGTACTTCTCCTTTTCCTTTTCATCGCTACCGTACTTCTCCTTTTCCTTTTCATCGTTACCGTTCTTCTCCTTTTCCTTTTCATCGATACCGTACTTCTCCTTTTCCTTTTCATCGTTACCGTACTTCTCCTTTTCCTTTTCATCGCTACCGTACTTCTCCTTTTCCTTTTCATCGCTACCGTACTTCTCCTTTTCCTTTTCATCGTTACCGTACTTCTCCTTTTCCTTTTCATCGCTACCGTACTTCTCCTTTTCCTTTTCATCGTTACCGTTCTTCTCCTTTTCCTTTTCATCGCTACCGTACTTCTCCTTTTCCTTTTCATCGTTACCGTACTTCTCCTTTTCCTTTTCATCGTTACCGTACTTCTCCTTTTCCTTTTCATCGTTACCGTTCTTCTTTTCCTTTTCATCACTACCGTACTTCTCCTTTTCCTTTTCATCGTTACCGTACTTctccttttcctttttatcGTTACCGTACTTTTCCTTTTCCTTTTCATCGTTACCGTACTTCTCCTTTTCCTTTTCATCGTTACCGTACTTCTCCTTTTCCTTTTCATCGTTACCGTACTTCTCCTTTTCCTTTTCATCGTTACCGTACTTCTCCTTTTCCTTTTCATCGCTACCGTACTTCTTCTCCTCTTCTTCTCATCGCTACCGTATTTCTTTGGCGGAGGGCACGTTGCCGGCTCTCCAGGG'.upper()
consensus = complement_dna(consensus_rev)[::-1]

start_rev = 'TGGCGGAGGGCACGTTGCCGGCTCTCCAGGG'
start = complement_dna(start_rev)[::-1]
end_rev = 'CAAACCTGCGCACAGCACGGAGACGATGCTC'
end = complement_dna(end_rev)[::-1]

with open('/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/truncated_reads.fasta', 'w') as out:
    # Use the semiglobal alignment to find the start and end in the reads

    for read, read_id in zip(tqdm(reads),read_ids):



        # Look at all the possible configurations of the read to find the ones where start and end fit the best
        seq_rev = read.upper()
        seq = complement_dna(seq_rev)[::-1]


        x, end_pos, end_score = find_in_seqs([seq], end)
        x, start_pos_rev, start_score = find_in_seqs([seq_rev], start_rev)


        # x, end_pos, x = find_in_seqs([seq], consensus)
        # x, start_pos_rev, x = find_in_seqs([seq_rev], consensus_rev)

        start_pos = len(seq) - start_pos_rev


        if end_pos < start_pos or end_score/len(end) < 0.5 or start_score/len(start) < 0.5:
            continue


        # seqs =  [read, read[::-1], complement_dna(read), complement_dna(read[::-1])]
        # start_mats = []
        # end_mats = []
        #
        # seq1, start_pos, start_score = find_in_seqs(seqs, consensus)
        # seq2, end_pos, end_score = find_in_seqs(seqs,end)

        # if not seq1 == seq2: break

        # Extract the region
        try:
            # HVD = seq1[start_pos - len(consensus) - 30: start_pos + 30]
            HVD = seq[start_pos :end_pos]
        except:
            print('read to short')


        # print everything

        out.write(read_id + '\n')
        out.write(HVD + '\n')
