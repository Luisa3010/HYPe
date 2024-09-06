##### Helper script to replace the sequence ids by their DNA sequences in Lukas files

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

infile = "/home/luisa/Downloads/Mappe13.csv"
outfile = "/home/luisa/Downloads/Mappe14.csv"
with open(infile, 'r') as file, open(outfile, 'w') as out:
    for i, line in enumerate(file):
        original_line = line.strip('\n')
        line = (',').join(line.split(',')[1:-1])
        for motif_aa, motif_dna in motifs.items():
            line = line.replace(motif_aa, motif_dna)
        out.write(original_line +',' + line + '\n')





