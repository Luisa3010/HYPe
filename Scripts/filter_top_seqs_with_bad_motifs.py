###### Look at the no amp reads and get those sequences that can not possibly be
# from the putative germline alleles



class Seq(object):

    motifs = []

    def __init__(self, motifs: list):
        self.motifs = motifs

    def __eq__(self, other):
        if not len(self.motifs) == len(other.motifs):
            return False

        for motif1, motif2 in zip(self.motifs, other.motifs):
            if not (motif1 == motif2 or motif1 == 'X' or motif2 == 'X'):
                return False

        return True


    def __str__(self):
        string = ''
        for motif in self.motifs:
            string += motif + ' '
        return string[:-1]



top_seqs = ['YERGGG SDRGGG RDRGD YERGGG SNRGGG SNRGGG SNRGGG SNRGGG RDRGD YERGGG SNRGGG SDRGD YERGGG SDRGGG RDRGD YERGGG SNRGGG SNRGGG REGGD YERGGG SNRGGG RDRGD YERGGG',
            'YERGGG SDRGGG RDNKRG SDRGD YERGGG SDRGGG SNRGGG SDRGD YERGGG SNRGGG SNRGGG SDRGD YERGGG SNRGGG SDRGD YERGGG SNRGGG SNRGGG SNRGGG REGGD YERGGG SNRGGG RDRGD YERGGG',
            'YERGGG RDNKRG SNRGGG RDRGD YERGGG SDRGGG SDRGGG RDRGD YERGGG SNRGGG SNRGGG RDRGD YERGGG SNRGGG SDRGD YERGGG SNRGGG SDRGD YERGGG SNRGGG RDRGGG REGGD YERGGG SNRGGG RDRGD YERGGG',
            'YERGGG SDRGGG RDRGD YERGGG SNRGGG SNRGGG RDRGD YERGGG SNRGGG SDRGD YERGGG SNRGGG SNRGGG SDRGD YERGGG SNRGGG SNRGGG REGGD YERGGG SNRGGG RDRGD YERGGG',
            'YERGGG SDRGGG RDNKRG SNRGGG SDRGD YERGGG SNRGGG REGGD YERGGG SNRGGG RDRGD YERGGG SNRGGG RDRGD YERGGG SNRGGG SDRGE YERGGG SNRGGG SDRGD YERGGG SNRGGG RDRGD YERGGG'
            'YERGGG SDRGGG RDNKRG SNRGGG RDRGD YERGGG SNRGGG SDRGD YERGGGG SNRGGG SDRGD YERGGG SNRGGG SDRGD YERGGG SNRGGG RDRGD YERGGG']


path = "/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_CAS9/Cas9_nanopore_round1_Newton_minimap2_clipped_tolerance5_sorted_HYP1_region_merged_q17_aa_bTrue_diff0"
seqs = []

# get all the reads
with open(path, 'r') as file:
    for line in file:
        if not line.startswith('>'):
            seqs.append(Seq(line.strip('\n').split(' ')))



for seq in seqs:
    is_top = False
    for top_seq in top_seqs:
        if Seq(top_seq.split(' ')) == seq:
            is_top = True
    if not is_top:
        print(seq)


