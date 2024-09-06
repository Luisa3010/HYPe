###### for a list of blocks, get those sequences that are not made from these
# blocks (or from SNPs) to see what they are made from


def check_only_SNP(string1, string2):
    if len(string1) == len(string2):
        count_diffs = 0
        for a, b in zip(string1, string2):
            if a!=b:
                if count_diffs > 0: return False
                count_diffs += 1
        return True
    else:
        return False


blocks = ["YERGGG SDRGGG RDRGD",
          "YERGGG SNRGGG SNRGGG SNRGGG SNRGGG RDRGD",
          "YERGGG SNRGGG SDRGD",
          "YERGGG SNRGGG SNRGGG REGGD",
          "YERGGG SNRGGG RDRGD",
          "YERGGG SDRGGG RDNKRG SDRGD",
          "YERGGG SDRGGG SNRGGG SDRGD",
          "YERGGG SNRGGG SNRGGG SDRGD",
          "YERGGG SNRGGG SNRGGG SNRGGG REGGD",
          "YERGGG RDNKRG SNRGGG RDRGD",
          "YERGGG SDRGGG SDRGGG RDRGD",
          "YERGGG SNRGGG SNRGGG RDRGD",
          "YERGGG SNRGGG RDRGGG REGGD",
          "YERGGG SNRGGG RDRGD"
          ]

filepath = "/home/luisa/Documents/Work/Uni_Cambridge/Data/working_files/tails1_noblock_sequences"

not_SNPs = []
SNPs = []
with open(filepath, 'r') as file:
    for line in file:
        is_SNP = False
        for block in blocks:
            is_SNP = is_SNP or check_only_SNP(line.strip("\n"),block)

        if not is_SNP:
            not_SNPs.append(line.strip("\n"))
        else:
            SNPs.append(line.strip("\n"))


print(f"SNPS: {len(SNPs)}")
for SNP in SNPs:
    print(SNP)


print(f"\n\nNon-SNPS: {len(not_SNPs)}")
for noSNP in not_SNPs:
    print(noSNP)