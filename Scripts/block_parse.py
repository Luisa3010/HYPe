###### Parse the HVDs of a file into blocks to explore if there can be a
# structure found in this


blocks = {
    'YERGGG SDRGGG RDRGD': 1,
    'YERGGG SNRGGG SNRGGG SNRGGG SNRGGG RDRGD': 2,
    'YERGGG SNRGGG SDRGD': 3,
    'YERGGG SNRGGG SNRGGG REGGD': 4,
    'YERGGG SNRGGG RDRGD': 5,
    'YERGGG\n': 6,
    'YERGGG SDRGGG RDNKRG SDRGD': 7,
    'YERGGG SDRGGG SNRGGG SDRGD': 8,
    'YERGGG SNRGGG SNRGGG SDRGD': 9,
    'YERGGG SNRGGG SNRGGG SNRGGG REGGD': 10,
    'YERGGG RDNKRG SNRGGG RDRGD': 11,
    'YERGGG SDRGGG SDRGGG RDRGD': 12,
    'YERGGG SNRGGG SNRGGG RDRGD': 13,
    'YERGGG SNRGGG RDRGGG REGGD': 14,
    'YERGGG SDRGGG RDNKRG SNRGGG RDRGD': 15,
    'YERGGG SDRGGG RDNKRG SNRGGG SDRGD': 16,
    'YERGGG SNRGGG REGGD': 17,
    'YERGGG SNRGGG SDRGE': 18
}
for i in range(1,9):
    infile = f'/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_semiglobal/singlecysts{i}.phred15.400to1400bp.fastq.sam.alignscore2000.fulllength_aa_diff0'

    n = 0

    with open(infile, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                line = line.replace('YQRGGG', 'YERGGG').replace('RDDQRG',
                                                                   'RDNKRG').replace(
                    'GNRGGG', 'SNRGGG').replace('SSRGD','SDRGD')

                for block,num in blocks.items():
                    line = line.replace(block,str(num))
                if line.replace(' ', '').isdigit():
                    print(line)
                    n +=1


    # print(n)