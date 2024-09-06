##### Get the counts for the top most frequent alleles in all reads from HYP1 G. Pallida


import os.path
from collections import Counter

top_seqs = [
        'YERGGG SDRGGG RDRGD YERGGG SNRGGG SNRGGG SNRGGG SNRGGG RDRGD YERGGG SNRGGG SDRGD YERGGG SDRGGG RDRGD YERGGG SNRGGG SNRGGG REGGD YERGGG SNRGGG RDRGD YERGGG',
        'YERGGG SDRGGG RDNKRG SDRGD YERGGG SDRGGG SNRGGG SDRGD YERGGG SNRGGG SNRGGG SDRGD YERGGG SNRGGG SDRGD YERGGG SNRGGG SNRGGG SNRGGG REGGD YERGGG SNRGGG RDRGD YERGGG',
        'YERGGG RDNKRG SNRGGG RDRGD YERGGG SDRGGG SDRGGG RDRGD YERGGG SNRGGG SNRGGG RDRGD YERGGG SNRGGG SDRGD YERGGG SNRGGG SDRGD YERGGG SNRGGG RDRGGG REGGD YERGGG SNRGGG RDRGD YERGGG',
        'YERGGG SDRGGG RDRGD YERGGG SNRGGG SNRGGG RDRGD YERGGG SNRGGG SDRGD YERGGG SNRGGG SNRGGG SDRGD YERGGG SNRGGG SNRGGG REGGD YERGGG SNRGGG RDRGD YERGGG',
        'YERGGG SDRGGG RDNKRG SNRGGG SDRGD YERGGG SNRGGG REGGD YERGGG SNRGGG RDRGD YERGGG SNRGGG RDRGD YERGGG SNRGGG SDRGE YERGGG SNRGGG SDRGD YERGGG SNRGGG RDRGD YERGGG',
        'YERGGG SDRGGG RDNKRG SNRGGG RDRGD YERGGG SNRGGG SDRGD YERGGG SNRGGG SDRGD YERGGG SNRGGG SDRGD YERGGG SNRGGG RDRGD YERGGG' # from singecyst3

        ]

directory = '/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_semiglobal/grouped'
outpath = '/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/Top_Counts_by_file_Heads_Tails.csv'


files = []
for file in os.listdir(directory):
    if os.path.isfile(os.path.join(directory, file)) and ('single' in file or 'head' in file or 'tails' in file) and 'aa' in file and not ('watch' in file or 'petri' in file or 'ref' in file or 'neg' in file or 'plasmid' in file ):
    # if os.path.isfile(os.path.join(directory, file)) and ('heads' in file or 'tails' in file) and 'aa' in file and not ('watch' in file or 'petri' in file or 'ref' in file or 'neg' in file ):
        files.append(file)


with open(outpath, 'w') as out:
    # out.write("File,Sequence,Count,Total,Total Top\n")
    for path in files:

        seqs = []
        n = 0
        top_n = 0
        top_unique_n = 0
        # get the sequences
        with open(os.path.join(directory, path), 'r') as file:
            for line in file:
                if not line.startswith('>'):
                    seqs.append(line.strip('\n'))
                    n += 1

            # Get the counts
            counts = Counter(seqs)
            # get the Counts of the top sequences
            # counts = [(s,counts[s]) for s in top_seqs]
            for seq, count in counts.items():
                top_n += count
                if count/n > 0.01:
                    top_unique_n += 1


            category = "Heads" if 'heads' in path else ("Tails" if "tails" in path else ('Singlecyst' if 'single' in path else ''))
            # print( "Heads " + str(top_unique_n) if 'heads' in path else "tails " + str(top_unique_n))
            print(category +' '  + str(top_unique_n))
            # for i, (seq, count) in enumerate(counts):
            #     out.write(f"{path},{seq},{count},{n},{top_n} \n")

