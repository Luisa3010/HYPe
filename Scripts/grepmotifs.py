##### grep but for a lot of strings and all the output are written into one file

import os
import argparse

def count_strings_in_files(text_files_folder, strings_file):
    # Read the strings to be counted from the strings file
    with open(strings_file, 'r') as sf:
        strings_to_count = [line.strip() for line in sf.readlines()]

    # Initialize a dictionary to keep the count of each string
    string_counts = {string: 0 for string in strings_to_count}

    # Iterate over all text files in the specified folder
    for filename in os.listdir(text_files_folder):

        with open(os.path.join(text_files_folder, filename), 'r') as tf:
            content = tf.read()
            for string in strings_to_count:
                string_counts[string] += content.count(string)

        print(f'{filename}:')
        for string, count in string_counts.items():
            print(f'{string}: {count}')
        string_counts = {string: 0 for string in strings_to_count}

    return string_counts

def main():
    parser = argparse.ArgumentParser(description="Count occurrences of specific strings in text files.")
    parser.add_argument('-i1', '--input_folder', required=True, help='Folder containing all the text files.')
    parser.add_argument('-i2', '--strings_file', required=True, help='File with strings to be counted.')
    args = parser.parse_args()

    string_counts = count_strings_in_files(args.input_folder, args.strings_file)

    # Print the results to stdout
    # for string, count in string_counts.items():
    #     print(f'{string}: {count}')

if __name__ == "__main__":
    main()
