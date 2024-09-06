#### helper file to create a sankey file to visualize reads - outdated, use visualize_sequences_as_sankey.py instead

import re


def insert_marker_and_reformat(input_file, output_file):
    with open(input_file, 'r') as file:
        content = file.read()

    # Insert [1] between numbers and digits
    pattern = re.compile(r'(\d+)(\S)')
    modified_content = pattern.sub(r'\1[1]\2', content)

    # Split the content by the marker to place each part on a new line
    parts = modified_content.split('[1]')

    # Remove trailing whitespace from each part and join them with newlines
    formatted_content = '\n'.join(
            [parts[i] + '[1]' + parts[i + 1] for i in range(len(parts) - 1)])

    with open(output_file, 'w') as file:
        file.write(formatted_content)


# Example usage
input_file = ('tmp')
output_file = 'output.txt'
insert_marker_and_reformat(input_file, output_file)
