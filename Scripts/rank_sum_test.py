##### perform a Wilcoxon rank sum test (man whitney u test) to get statistical
# significant differences between two non parametric distributed groups


import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency
import argparse
from scipy.stats import ranksums
from scipy.stats import ks_2samp




def main():
    # Combine data into a DataFrame

    # get the arguments from the command line
    parser = argparse.ArgumentParser()
    parser.add_argument('-i1', '--input1', required = True, help = 'List of numbers from category A')
    parser.add_argument('-i2', '--input2', required = True, help = 'List of numbers from category B')
    parser.add_argument('-c1', '--input3', required = False, default = '', help = 'Counts of category A')
    parser.add_argument('-c2', '--input4', required = False, default = '', help = 'Counts of category B')

    args = parser.parse_args()
    infile1 = args.input1
    infile2 = args.input2
    infile3 = args.input3
    infile4 = args.input4

    category_a = []
    with open(infile1, 'r') as file:
        for line in file:
            category_a.append(float(line.strip()))
            # category_a.append(line.strip()) # uncomment for text data


    category_b =[]
    with open(infile2, 'r') as file:
        for line in file:
            category_b.append(float(line.strip()))
            # category_b.append(line.strip()) # uncomment for text data

    counts_a = []
    if not infile3 == '': 
        with open(infile3, 'r') as file:
            for i, line in enumerate(file):
                for j in range(int(line.strip())):
                    counts_a.append(category_a[i]) 
        category_a = counts_a

    counts_b = []
    if not infile4 == '':
        with open(infile4, 'r') as file:
            for i, line in enumerate(file):
                for j in range(int(line.strip())):
                    counts_b.append(category_b[i])
        category_b = counts_b


            # Perform the Wilcoxon rank-sum test
    stat, p = ranksums(category_a, category_b)

    # Output the results
    print(f"Statistic: {stat}")
    print(f"P-Value: {p}")
    if p < 0.05:
        print(
            "The distributions of categories A and B are significantly different.")
    else:
        print(
            "There is no significant difference \nbetween the distributions of categories A and B.")


if __name__ == '__main__':
    main()