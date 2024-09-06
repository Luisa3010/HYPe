###### check if barcodes for certain sequences containing only YERGGG as HVD
# and which seem to be in the wrong file are correct

ids = [272,
       1028,
       2232,
       3436,
       3444,
       4032,
       4264,
       4478,
       4784,
       4878,
       5052,
       5088,
       5110,
       5276,
       5384,
       5680,
       5890,
       6036,
       6428,
       6954,
       8166,
       8612,
       8990,
       9148,
       9780,
       9916,
       10974,
       11784,
       11828,
       11978,
       12660,
       13378,
       14150]


barcodes_fwd = [
    "CAGGTTACTCCTCCGTGAGTCTGA",
    "AACCAAGACTCGCTGTGCCTAGTT",
    "TGAGAGACAAGATTGTTCGTGGAC",
    "GTTCCTCGTGCAGTGTCAAGAGAT",
    "ACCACAGGAGGACGATACAGAGAA",
    "CGGATGAACATAGGATAGCGATTC",
    "TACAGTCCGAGCCTCATGTGATCT",
    "TGCGTACAGCAATCAGTTACATTG",
    "CATGTTCAACCAAGGCTTCTATGG",
    "AGAACGACTTCCATACTCGTGTGA"
    ]


barcodes_rev = [
    "TCAGACTCACGGAGGAGTAACCTG",
    "AACTAGGCACAGCGAGTCTTGGTT",
    "GTCCACGAACAATCTTGTCTCTCA",
    "ATCTCTTGACACTGCACGAGGAAC",
    "TTCTCTGTATCGTCCTCCTGTGGT",
    "GAATCGCTATCCTATGTTCATCCG",
    "AGATCACATGAGGCTCGGACTGTA",
    "CAATGTAACTGATTGCTGTACGCA",
    "CCATAGAAGCCTTGGTTGAACATG",
    "TCACACGAGTATGGAAGTCGTTCT"
]

rep1_mix1_bc_ids = [
        (1,10),
        (9,1),
        (9,2),
        (9,3),
        (9,4),
        (9,5),
        (9,6),
        (9,7),
        (9,8),
        (4,2),
        (4,3),
        (4,4),
        (4,5),
        (4,6),
        (4,7),
        (4,8),
        (4,9),
        (4,10),
        (5,2),
        (5,3),
        (5,4),
        (5,5),
        (5,6),
        (5,7),
        (5,8),
        (5,9),
        (5,10)
]
inputfile = "/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT/raw_data/PAS73765_pass_af667136_d703a591_0_sub1000.fasta"

seqs = {}

with open(inputfile, 'r') as file:
    prev = ''
    for i, line in enumerate(file):
        if i+1 in ids:
            print(prev[:-1])
            print(line)
            seqs[prev[:-1]] = line[:-1]
        prev = line


is_correct_mix = []
for id, seq in seqs.items():
    fwd, rev = ('', '')
    for i, barcode_fwd in enumerate(barcodes_fwd):
        if barcode_fwd in seq:
            fwd = i + 1
    for i, barcode_rev in enumerate(barcodes_rev):
        if barcode_rev in seq:
            rev = i + 1

    print(f'{id} fwd = {fwd}, rev = {rev}')
    if not fwd == '' and not rev == '':
        print((fwd,rev) in rep1_mix1_bc_ids)
        is_correct_mix.append((fwd,rev) in rep1_mix1_bc_ids)



print(f'Ratio of sequences with single YERGGG in correct plasmid solutions: {sum(is_correct_mix)/len(is_correct_mix)}')





