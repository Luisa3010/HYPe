######## Create a tree like diagram for a File of G. Pallida HYP1 HVDs outdated,
# use vsualize_sequences_as_sankey.py instead


import os.path
from tqdm import tqdm

import plotly.graph_objects as go

motif_colors= {
    'YERGGG': '#8a2424' ,
    'SNRGGG': '#00450c',
    'SDRGD': '#4a3200',
    'RDRGD': '#996800',
    'SDRGGG': '#00ba37',
    'RDNKRG': '#0a4b78',
    'REGGD': '#bd8600',
    'RDRGGG': '#211600',
    'SDRGE': '#dba617',
    'GNRGGG': '#003008',
    'RDDQRG': '#135e96',
    'END': '#ffffff'
}
# seqs = ['YERGGG SDRGGG RDRGD YERGGG SNRGGG SNRGGG SNRGGG SNRGGG RDRGD YERGGG SNRGGG SDRGD YERGGG SDRGGG RDRGD YERGGG SNRGGG SNRGGG REGGD YERGGG SNRGGG RDRGD YERGGG',
#         'YERGGG SDRGGG RDNKRG SDRGD YERGGG SDRGGG SNRGGG SDRGD YERGGG SNRGGG SNRGGG SDRGD YERGGG SNRGGG SDRGD YERGGG SNRGGG SNRGGG SNRGGG REGGD YERGGG SNRGGG RDRGD YERGGG',
#         'YERGGG SDRGGG RDNKRG SDRGD YERGGG SDRGGG SNRGGG SDRGD YERGGG SNRGGG SNRGGG SDRGD YERGGG SNRGGG SDRGD YERGGG SNRGGG SNRGGG SNRGGG REGGD YERGGG SNRGGG RDRGD YERGGG']

for l in tqdm(range(1,11)):
    infile = f'/home/luisa/Documents/Work/Uni_Cambridge/Data/HYP1_ONT_semiglobal/grouped/sorted/tails{l}.phred15.400to1400bp.alignscore2000.fulllength_aa_diff0'

    labels_list = []
    colors_list = []
    xs_list = []
    values_list = []
    max_id = 0


    with open(infile, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                seq = line.strip('\n')

                # get the motifs as labels
                labels = seq.split(' ')
                labels.append('END')
                # get the colors for each motif
                colors = [motif_colors[motif] for motif in labels]

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

    for i in range(len(labels_list)):
        j = 0
        while j < i:
            labels1 = labels_list[i] # current seq
            labels2 = labels_list[j] # seq for comparison
            # if they are exactly the same remove labels1 and make the values bigger by 1 for labels2
            if labels1 == labels2:
                for k in range(min(len(labels1), len(labels2))):
                    xs_list[i][k] = xs_list[j][k]
                break
            # compare each position, if they are the same, change the corresponding x for label1 to the x of label2
            for k in range(min(len(labels1),len(labels2))):
                if labels1[k] == labels2[k]:
                    xs_list[i][k] = xs_list[j][k]
                # if they are not the same move to the next seq
                else:
                    break
            j+=1


    labels = []
    colors = []
    sources = []
    targets = []
    values = []
    edge_colors = []

    for i in range(len(labels_list)):
            labels.extend(labels_list[i])
            colors.extend(colors_list[i])
            sources.extend(xs_list[i][:-1])
            targets.extend(xs_list[i][1:])
            values.extend(values_list[i])
            for j in range(len(xs_list[i])-2):
                edge_colors.append('lightgrey')
            edge_colors.append('rgba(1,1,1,0)')





    fig = go.Figure(data=[go.Sankey(
        textfont=dict(color="rgba(0,0,0,0)", size=1),
        node = dict(
          pad = 15,
          thickness = 20,
          line = dict(color = 'white', width = 0),
          label = labels,
          color = colors

        ),
        link = dict(
            source = sources,
            target = targets,
            value = values,
            color = edge_colors
      ))])


    title = f'Tails{l} - HVD Allele Variants'
    fig.update_layout(title_text=title, font_size=30)
    fig.write_image(f'/home/luisa/Documents/Work/Uni_Cambridge/Diagrams/{title}_sankey.png',scale=10, width=1680, height=1080)


    # Get legend of colors colors
    # make into neat little script that you can plug sequences in