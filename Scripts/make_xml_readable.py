##### take an xml file and convert it to be humanly readable

import xml.dom.minidom


def format_xml(input_file, output_file):
    with open(input_file, 'r') as f:
        xml_data = f.read()

    dom = xml.dom.minidom.parseString(xml_data)  # parse the string
    pretty_xml_as_string = dom.toprettyxml()  # format the XML string

    with open(output_file, 'w') as f:
        f.write(pretty_xml_as_string)


# Usage
path = "/home/luisa/Documents/Work/Uni_Cambridge/Data/plasmids/HYP1_plasmids_for_vincent/883334512  Dataset for C:/Users/Unnati/OneDrive - University of Cambridge/Data/Sequencing"
input_file = '/HYP1_TA_clones_all_exon2_PCR5_plus_KODX_clones'
output_file = '/formatted_output.xml'
format_xml(path + input_file, path + output_file)