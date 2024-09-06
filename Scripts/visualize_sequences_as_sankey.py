###### make a tree like visualization for HVDs for G. Pallida HYP1 that can be input through a user
# interface


import plotly.graph_objects as go
import tkinter as tk
from tkinter import filedialog

def open_file():
    file_path = filedialog.askopenfilename(filetypes=[("Text files", "*.txt"), ("All files", "*.*")])
    if file_path:
        with open(file_path, 'r') as file:
            content = file.read()
            input_text.delete("1.0", tk.END)
            input_text.insert(tk.END, content)

def create_viz():
    title = title_entry.get().strip()
    content = input_text.get("1.0", tk.END)

    motif_colors= {
        'END': '#ffffff',
        'END1': '#ffffff',
        # YERGGGs are red, YQRGGGs are violett
        'YERGGG1': '#8a2424',
        'YERGGG2': '#b32d2e',
        'YERGGG3': '#d63638',
        'YERGGG4': '#e65054',
        'YERGGG5': '#f86368',
        'YQRGGG1': '#8a245e',
        'YQRGGG2': '#b32d64',
        # SNRGGGS are moss green
        'SNRGGG1': '#00450c',
        'SNRGGG2': '#005c12',
        'SNRGGG3': '#007017',
        'SNRGGG4': '#008a20',
        'SNRGGG5': '#00a32a',
        # SDRGD is yellow
        'SDRGD1': '#996800',
        'SDRGD2': '#bd8600',
        'SSRGD1': '#dba617',
        # RDRGD is brown
        'RDRGD1': '#4a3200',
        # SDRGGG is light green
        'SDRGGG1': '#1ed14b',
        'SDRGGG2': '#68de7c',
        'SDRGGG3': '#b8e6bf',
        # RDNKRG is blue
        'RDNKRG1': '#0a4b78',
        # REGGD orange yellow
        'REGGD1': '#994000',
        # RDRGGG is grey
        'RDRGGG1': '#787c82',
        # SDRGE is yellow
        'SDRGE1': '#755100',
        # GNRGGGS are dark green
        'GNRGGG1': '#003008',
        'GNRGGG2': '#003008',
        'GNRGGG3': '#003008',
        # RDDQRg is lighter blue
        'RDDQRG1': '#135e96',
        # unknown motifs are black
        'X1': 'black'

    }

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
            "END": "END"

    }

    # Remove headers if the input is a fasta file
    content = content.strip('\n')
    seqs = [seq for seq in content.split('\n') if not seq.startswith('>')]

    # Filter for sequences by criteria
    # Remove duplicates
    # seqs = list(set(seqs))
    # top_seqs = ['YERGGG SDRGGG RDRGD YERGGG SNRGGG SNRGGG SNRGGG SNRGGG RDRGD YERGGG SNRGGG SDRGD YERGGG SDRGGG RDRGD YERGGG SNRGGG SNRGGG REGGD YERGGG SNRGGG RDRGD YERGGG',
    #             'YERGGG SDRGGG RDNKRG SDRGD YERGGG SDRGGG SNRGGG SDRGD YERGGG SNRGGG SNRGGG SDRGD YERGGG SNRGGG SDRGD YERGGG SNRGGG SNRGGG SNRGGG REGGD YERGGG SNRGGG RDRGD YERGGG',
    #             'YERGGG RDNKRG SNRGGG RDRGD YERGGG SDRGGG SDRGGG RDRGD YERGGG SNRGGG SNRGGG RDRGD YERGGG SNRGGG SDRGD YERGGG SNRGGG SDRGD YERGGG SNRGGG RDRGGG REGGD YERGGG SNRGGG RDRGD YERGGG',
    #             'YERGGG SDRGGG RDRGD YERGGG SNRGGG SNRGGG RDRGD YERGGG SNRGGG SDRGD YERGGG SNRGGG SNRGGG SDRGD YERGGG SNRGGG SNRGGG REGGD YERGGG SNRGGG RDRGD YERGGG',
    #             'YERGGG SDRGGG RDNKRG SNRGGG SDRGD YERGGG SNRGGG REGGD YERGGG SNRGGG RDRGD YERGGG SNRGGG RDRGD YERGGG SNRGGG SDRGE YERGGG SNRGGG SDRGD YERGGG SNRGGG RDRGD YERGGG'
    #             ]
    # seqs = [seq for seq in seqs if len(seq.split(' ')) > 2 and seq.split(' ')[2] == 'RDRGD']
    # seqs = [seq for seq in seqs if not seq in top_seqs]
    # seqs = [seq for seq in seqs if len(seq.split(' ')) < 21 ]
    # seqs = top_seqs


    labels_list = []
    colors_list = []
    xs_list = []
    values_list = []
    max_id = 0


    for seq in seqs:

        # get the motifs as labels
        labels = seq.split(' ')
        labels.append('END')
        # Determine which version of the motifs names is present and get the colors for each motif
        if labels[0] in motif_colors.keys():
            colors = [motif_colors[motif] for motif in labels]
        elif labels[0] in motifs.keys():
            colors = [motif_colors[motifs[motif]] for motif in labels]
        elif labels[0] + '1' in motif_colors.keys():
            colors = [motif_colors[motif + '1'] for motif in labels]


        # append the indices of the labels
        x = []
        value = []
        for i in range(max_id, max_id + len(labels)):
            x.append(i)
            value.append(1)
        max_id = i + 1

        labels_list.append(labels)
        colors_list.append(colors)
        xs_list.append(x)
        values_list.append(value)

    # Check for the overlap of sequences
    for i in range(len(labels_list)):
        j = 0
        while j < i:
            labels1 = labels_list[i] # current seq
            labels2 = labels_list[j] # seq for comparison
            # if they are exactly the same duplicate the indices from the first entry
            if labels1 == labels2:
                for k in range(min(len(labels1), len(labels2))):
                    xs_list[i][k] = xs_list[j][k]
                break
            # Compare each position, if they are the same, change the corresponding x for label1 to the x of label2

            ##### Comment this section if you want a non-treelike structure #####
            for k in range(min(len(labels1),len(labels2))):
                if labels1[k] == labels2[k]:
                    xs_list[i][k] = xs_list[j][k]
                # If they are not the same move to the next seq
                elif not (labels1[k] == 'X' or labels2[k] == 'X') :
                    break
                # else:
                #     break
            ##### End of section to comment out #####

            j+=1



    labels = []
    colors = []
    sources = []
    targets = []
    values = []
    edge_colors = []

    # Convert into format required to create the Sankey diagram
    for i in range(len(labels_list)):
            labels.extend(labels_list[i])
            colors.extend(colors_list[i])
            sources.extend(xs_list[i][:-1])
            targets.extend(xs_list[i][1:])
            values.extend(values_list[i])

            # Append grey edges between the motifs
            for j in range(len(xs_list[i])-2):
                edge_colors.append('lightgrey')

            # Append transparent edge for the end
            edge_colors.append('rgba(1,1,1,0)')

    # Create the sankey
    fig = go.Figure(data=[go.Sankey(
        textfont=dict(color="rgba(0,0,0,0)", size=1),
        node = dict(
          pad = 20,
          thickness =30,
          line = dict( width = 0),
          label = labels,
          color = colors

        ),
        link = dict(
            source = sources,
            target = targets,
            value = values,
            color = edge_colors
      ))])


    fig.update_layout(title_text=title, font_size=30)
    # fig.write_image(f'/home/luisa/Documents/Work/Uni_Cambridge/Diagrams/{title}_sankey.png',scale=10, width=1680, height=1080)
    fig.show()


# Create the main window
root = tk.Tk()
root.title("Input Sequences")

# Add a label and entry for the title
title_label = tk.Label(root, text="Input Title for your diagram:")
title_label.pack(pady=5)

title_entry = tk.Entry(root, width=50)
title_entry.pack(pady=5)


# Add a label with text above the input box
label = tk.Label(root, text="Please enter text or open a file:")
label.pack(pady=5)

# Create a Text widget
input_text = tk.Text(root, height=10, width=50)
input_text.pack(pady=10)

# Create a Button to open a file
file_button = tk.Button(root, text="Open File", command=open_file)
file_button.pack(pady=5)

# Create a Button to get the input text
input_button = tk.Button(root, text="Create Viz", command=create_viz)
input_button.pack(pady=5)

# Run the application
root.mainloop()

