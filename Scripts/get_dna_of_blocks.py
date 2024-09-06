######## Look at what blocks are present in the data and at which frequency they occur

motifs = {
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

# open file
path = '/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_semiglobal'

blocks_count = {}

for i in range(1,11):
    infile = path + f"/heads{i}.phred15.400to1400bp.alignscore2000.fulllength_dna_diff0"


    # read line
    with open(infile,'r') as file:
        for line in file:
            if not line.startswith('>'):
                # split at YERGGG/YQRGGG
                blocks = [block for block in line.strip('\n').split(' TAT')]
                # count the frequency of the block
                for block in blocks:
                    if not block.startswith('TAT'):
                        block = 'TAT' + block
                    # aa_block = (' ').join([motifs[dna_motif] for dna_motif in block.split(' ')])
                    aa_block = block
                    if aa_block in blocks_count.keys():
                        blocks_count[aa_block] += 1
                    else:
                        blocks_count[aa_block] = 1




# print the blocks with their counts
print(f"Block, Count")
for block, count in blocks_count.items():
    # if count > 100:
    print(f"{block}, {count}")
    # print(block)