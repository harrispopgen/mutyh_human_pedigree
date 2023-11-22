import os
import numpy as np
import pandas as pd
import zipfile
import re
from collections import Counter

from pybedtools import BedTool

import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import matplotlib.patches as patches

import dash
from jupyter_dash import JupyterDash
from dash import dcc, html
from dash.dependencies import Input, Output

import plotly.graph_objects as go
from plotly.offline import plot
import plotly.express as px

def sort_chromosomes(chrom):
    if chrom.startswith('chr'):
        chrom = chrom[3:]
    try:
        return int(chrom)
    except ValueError:
        return chrom
    
def parse_bed_file(filename):
    with open(filename, 'r') as file:
        intervals = []
        for line in file:
            parts = line.strip().split()
            intervals.append((parts[0], int(parts[1]), int(parts[2])))  # chrom, start, end
    return intervals

def intersect_bed_files(comparison_groups, bed_directory='dense_regions_beds'):
    """Intersect the BED files for each comparison group and save the result."""
    for group in comparison_groups:
        # Initialize with the first BED file
        base_bed = BedTool(os.path.join(bed_directory, f"{group[0]}_dense_clusters_sorted.bed"))
        
        # Intersect the base with subsequent files in the group
        for file_title in group[1:]:
            file_path = os.path.join(bed_directory, f"{file_title}_dense_clusters_sorted.bed")
            base_bed = base_bed.intersect(BedTool(file_path))

        # Save the intersected result as a new BED file
        output_filename = "_vs_".join(group) + "_overlap.bed"
        base_bed.saveas(os.path.join(bed_directory, output_filename))

C11askid = pd.read_csv("C11askid.csv")
C12askid = pd.read_csv("C12askid.csv")

C31askid = pd.read_csv("C31askid.csv")
C32askid = pd.read_csv("C32askid.csv")

C22asfather_C21askid = pd.read_csv("C22asfather_C21askid.csv")
C23asfather_C21askid = pd.read_csv("C23asfather_C21askid.csv")
C21asfather_C22askid = pd.read_csv("C21asfather_C22askid.csv")
C23asfather_C22askid = pd.read_csv("C23asfather_C22askid.csv")
C21asfather_C23askid = pd.read_csv("C21asfather_C23askid.csv")
C22asfather_C23askid = pd.read_csv("C22asfather_C23askid.csv")

C42asfather_C41askid = pd.read_csv("C42asfather_C41askid.csv")
C41asfather_C42askid = pd.read_csv("C41asfather_C42askid.csv")

P2P3asparents_P1askid = pd.read_csv("P2P3asparents_P1askid.csv")
P2P4asparents_P1askid = pd.read_csv("P2P4asparents_P1askid.csv")
P3P4asparents_P1askid = pd.read_csv("P3P4asparents_P1askid.csv")

P1P3asparents_P2askid = pd.read_csv("P1P3asparents_P2askid.csv")
P1P4asparents_P2askid = pd.read_csv("P1P4asparents_P2askid.csv")
P3P4asparents_P2askid = pd.read_csv("P3P4asparents_P2askid.csv")

P1P4asparents_P3askid = pd.read_csv("P1P4asparents_P3askid.csv")
P1P2asparents_P3askid = pd.read_csv("P1P2asparents_P3askid.csv")
P2P4asparents_P3askid = pd.read_csv("P2P4asparents_P3askid.csv")

P1P2asparents_P4askid = pd.read_csv("P1P2asparents_P4askid.csv")
P1P3asparents_P4askid = pd.read_csv("P1P3asparents_P4askid.csv")
P2P3asparents_P4askid = pd.read_csv("P2P3asparents_P4askid.csv")

# List of data frames and their titles
data_frames = [
    (C11askid, "C11_as_kid"),
    (C12askid, "C12_as_kid"),
    (C31askid, "C31_as_kid"),
    (C32askid, "C32_as_kid"),
    (C22asfather_C21askid, "C22asfather_C21askid"),
    (C23asfather_C21askid, "C23asfather_C21askid"),
    (C21asfather_C22askid, "C21asfather_C22askid"),
    (C23asfather_C22askid, "C23asfather_C22askid"),
    (C21asfather_C23askid, "C21asfather_C23askid"),
    (C22asfather_C23askid, "C22asfather_C23askid"),
    (C42asfather_C41askid, "C42asfather_C41askid"),
    (C41asfather_C42askid, "C41asfather_C42askid"),
    (P2P3asparents_P1askid, "P2P3asparents_P1askid"),
    (P2P4asparents_P1askid, "P2P4asparents_P1askid"),
    (P3P4asparents_P1askid, "P3P4asparents_P1askid"),
    (P1P3asparents_P2askid, "P1P3asparents_P2askid"),
    (P1P4asparents_P2askid, "P1P4asparents_P2askid"),
    (P3P4asparents_P2askid, "P3P4asparents_P2askid"),
    (P1P4asparents_P3askid, "P1P4asparents_P3askid"),
    (P1P2asparents_P3askid, "P1P2asparents_P3askid"),
    (P2P4asparents_P3askid, "P2P4asparents_P3askid"),
    (P1P2asparents_P4askid, "P1P2asparents_P4askid"),
    (P1P3asparents_P4askid, "P1P3asparents_P4askid"),
    (P2P3asparents_P4askid, "P2P3asparents_P4askid")
]

comparison_groups = [
    ['C22asfather_C21askid', 'C23asfather_C21askid'],
    ['C21asfather_C22askid', 'C23asfather_C22askid'],
    ['C21asfather_C23askid', 'C22asfather_C23askid'],
    ['P2P3asparents_P1askid', 'P2P4asparents_P1askid', 'P3P4asparents_P1askid'],
    ['P1P3asparents_P2askid', 'P1P4asparents_P2askid', 'P3P4asparents_P2askid'],
    ['P1P4asparents_P3askid', 'P1P2asparents_P3askid', 'P2P4asparents_P3askid'],
    ['P1P2asparents_P4askid', 'P1P3asparents_P4askid', 'P2P3asparents_P4askid']
]

comparison_titles = [
    'C21_mutations',
    'C22_mutations',
    'C23_mutations',
    'P1_mutations',
    'P2_mutations',
    'P3_mutations',
    'P4_mutations'
]

chr_lengths = {
    'chr1': 248956422,
    'chr2': 242193529,
    'chr3': 198295559,
    'chr4': 190214555,
    'chr5': 181538259,
    'chr6': 170805979,
    'chr7': 159345973,
    'chr8': 145138636,
    'chr9': 138394717,
    'chr10': 133797422,
    'chr11': 135086622,
    'chr12': 133275309,
    'chr13': 114364328,
    'chr14': 107043718,
    'chr15': 101991189,
    'chr16': 90338345,
    'chr17': 83257441,
    'chr18': 80373285,
    'chr19': 58617616,
    'chr20': 64444167,
    'chr21': 46709983,
    'chr22': 50818468
}

windowSize = 15000000
stepSize = 3000000
threshold = 4

# Iterate through data frames and their titles
for index, (mutation_data, title) in enumerate(data_frames):
    mutation_data['ref_alt_value'] = mutation_data['REF.ALT']
    mutation_data['Window Mutation Count'] = 0
    
    # Iterate through each chromosome
    for chrom, chr_data in mutation_data.groupby('X.CHROM'):
        chrLen = chr_lengths[chrom]  # Get the length of the current chromosome

        windowStart = 0
        windowEnd = windowSize

        while windowStart < chrLen:
            if windowEnd > chrLen:
                windowEnd = chrLen
            
            # Current window
            curr_window_data = chr_data[(chr_data['POS'] >= windowStart) & (chr_data['POS'] < windowEnd)]
            curr_window_count = len(curr_window_data)

            # Previous window (if it exists)
            prev_window_start = windowStart - stepSize
            prev_window_data = chr_data[(chr_data['POS'] >= prev_window_start) & (chr_data['POS'] < windowStart)]
            prev_window_count = len(prev_window_data) if prev_window_start >= 0 else 0

            # Next window
            next_window_end = windowEnd + stepSize
            next_window_data = chr_data[(chr_data['POS'] >= windowEnd) & (chr_data['POS'] < next_window_end)]
            next_window_count = len(next_window_data)

            if curr_window_count > threshold:
                for idx in curr_window_data.index:
                    mutation_data.at[idx, 'Window Mutation Count'] = curr_window_count

            # Move the window
            windowStart += stepSize
            windowEnd += stepSize

    mutation_data['Mutation Type'] = np.where(mutation_data['Window Mutation Count'] > threshold, 'Dense', 'Sparse')
    mutation_data['Cluster Distance'] = mutation_data.groupby('X.CHROM')['POS'].diff().fillna(0)
    
    mutation_data['Cluster ID'] = (mutation_data['Mutation Type'] == 'Dense') & (mutation_data['Cluster Distance'] > 3000000)
    mutation_data['Cluster ID'] = mutation_data.groupby('X.CHROM')['Cluster ID'].cumsum()

    clusters_grouped = mutation_data[mutation_data['Mutation Type'] == 'Dense'].groupby(['X.CHROM', 'Cluster ID'])
    cluster_positions = clusters_grouped.agg({'POS': ['min', 'max']})
    cluster_positions.columns = ['Start', 'End']
    
    # Create BED file
    bed_lines = []
    for (chrom, _), cluster_info in cluster_positions.iterrows():
        start = int(cluster_info['Start'])
        end = int(cluster_info['End'])
        bed_line = f"{chrom}\t{start}\t{end}\n"
        bed_lines.append(bed_line)

    bed_lines.sort(key=lambda line: (sort_chromosomes(line.split('\t')[0]), int(line.split('\t')[1])))
        
    # Save the BED file to the "dense_regions_beds" directory
    bed_file_name = f'{title}_dense_clusters_sorted.bed'
    bed_file_path = os.path.join('dense_regions_beds', bed_file_name)
    with open(bed_file_path, 'w') as bed_file:
        bed_file.writelines(bed_lines)
    
    
    filtered_mutations = mutation_data[mutation_data['Mutation Type'] == 'Sparse']
    
    # Save the CSV file to the "sparsity_filtered_mutations" directory
    csv_file_name = f'{title.replace(" ", "_")}_filtered_mutations.csv'
    csv_file_path = os.path.join('sparsity_filtered_mutations', csv_file_name)
    filtered_mutations.to_csv(csv_file_path, index=False)
    

intersect_bed_files(comparison_groups)

app = dash.Dash(__name__)

data_frame_selector = dcc.Dropdown(
    id='data-frame-selector',
    options=[{'label': title, 'value': i} for i, (_, title) in enumerate(data_frames)],
    value=0 
)

app.layout = html.Div([
    data_frame_selector,
    dcc.Graph(id='score-plot')
])

@app.callback(
    Output('score-plot', 'figure'),
    [Input('data-frame-selector', 'value')]
)

def update_plot(selected_index):
    mutation_data, title = data_frames[selected_index]

    clusters_grouped = mutation_data[mutation_data['Mutation Type'] == 'Dense'].groupby(['X.CHROM', 'Cluster ID'])
    cluster_positions = clusters_grouped.agg({'POS': ['min', 'max']})
    cluster_positions.columns = ['Start', 'End']

    chromosome_order = sorted(mutation_data['X.CHROM'].unique(), key=sort_chromosomes)

    start_traces = go.Scatter(
        x=cluster_positions['Start'],
        y=cluster_positions.index.get_level_values('X.CHROM'),
        mode='markers',
        marker=dict(
            size=10,
            symbol='triangle-up',
            color='green', 
        ),
        name='Start of Cluster',
    )

    end_traces = go.Scatter(
        x=cluster_positions['End'],
        y=cluster_positions.index.get_level_values('X.CHROM'),
        mode='markers',
        marker=dict(
            size=10,
            symbol='triangle-down',
            color='orange', 
        ),
        name='End of Cluster',
    )

    figure = {
        'data': [
            go.Scatter(
                x=mutation_data[mutation_data['Mutation Type'] == 'Dense']['POS'],
                y=mutation_data[mutation_data['Mutation Type'] == 'Dense']['X.CHROM'],
                mode='markers',
                marker={
                    'size': 10,
                    'opacity': 1.0,
                    'symbol': 'circle',  
                    'color': 'blue', 
                },
                name='Dense',
            ),
            go.Scatter(
                x=mutation_data[mutation_data['Mutation Type'] == 'Sparse']['POS'],
                y=mutation_data[mutation_data['Mutation Type'] == 'Sparse']['X.CHROM'],
                mode='markers',
                marker={
                    'size': 10,
                    'opacity': 1.0,
                    'symbol': 'circle-open',  
                    'color': 'red', 
                },
                name='Sparse',
            ),
            start_traces,
            end_traces,
        ],
        'layout': {
            'title': title + ' - Mutation Density Scores Across Chromosomes',
            'xaxis': {'title': 'Chromosome Position'},
            'yaxis': {'title': 'Chromosome', 'categoryorder': 'array', 'categoryarray': chromosome_order}
        }
    }

    return figure


if __name__ == '__main__':
    app.run_server(host='0.0.0.0', port=8050, debug=True)

# Folder to save HTML files
output_folder = "density_filtered_html_files"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

for idx, (mutation_data, title) in enumerate(data_frames):
    fig = update_plot(idx)  # You may need to adjust this part based on how you generate your figures
    plot(fig, filename=os.path.join(output_folder, f"{title.replace(' ', '_')}.html"), auto_open=False)

zip_file_name = "density_filtered_html_files.zip"

with zipfile.ZipFile(zip_file_name, 'w') as zipf:
    for root, _, files in os.walk(output_folder):
        for file in files:
            zipf.write(os.path.join(root, file), os.path.relpath(os.path.join(root, file), output_folder))

parent_mapping = {
    "C22asfather_C21askid": "C22",
    "C23asfather_C21askid": "C23",
    "C21asfather_C22askid": "C21",
    "C23asfather_C22askid": "C23",
    "C21asfather_C23askid": "C21",
    "C22asfather_C23askid": "C22",
    "P2P3asparents_P1askid": "P2P3",
    "P2P4asparents_P1askid": "P2P4",
    "P3P4asparents_P1askid": "P3P4",
    "P1P3asparents_P2askid": "P1P3",
    "P1P4asparents_P2askid": "P1P4",
    "P3P4asparents_P2askid": "P3P4",
    "P1P4asparents_P3askid": "P1P4",
    "P1P2asparents_P3askid": "P1P2",
    "P2P4asparents_P3askid": "P2P4",
    "P1P2asparents_P4askid": "P1P2",
    "P1P3asparents_P4askid": "P1P3",
    "P2P3asparents_P4askid": "P2P3"
}

def sort_comparison_chromosomes(chrom):
    # Split chromosome string into chromosome number and suffix (if present)
    parts = chrom.split('.')
    if parts[0].startswith('chr'):
        parts[0] = parts[0][3:]
    try:
        primary_key = int(parts[0])
    except ValueError:
        primary_key = parts[0]
    
    secondary_key = parts[1] if len(parts) > 1 else ''
    
    return (primary_key, secondary_key)

def update_chromosome_names(data, dataset_name):
    data = data.copy()  
    parent_identifier = parent_mapping.get(dataset_name, dataset_name)  
    data['X.CHROM'] = data['X.CHROM'].astype(str) + '.' + parent_identifier
    return data

comparison_group_selector = dcc.Dropdown(
    id='comparison_group_selector',
    options=[{'label': title, 'value': i} for i, title in enumerate(comparison_titles)],
    value=0 
)

app = dash.Dash(__name__)

app.layout = html.Div([
    comparison_group_selector,
    dcc.Graph(id='score-plot')
])

@app.callback(
    Output('score-plot', 'figure'),
    [Input('comparison_group_selector', 'value')]
)

def update_comparison_plot(selected_index):
    selected_comparison = comparison_groups[selected_index]
    title = comparison_titles[selected_index]
    combined_data = pd.concat([
        update_chromosome_names(data_frames[dataset_name_idx][0], dataset_name)
        for dataset_name in selected_comparison
        for dataset_name_idx, (_, dataset_name_title) in enumerate(data_frames)
        if dataset_name_title == dataset_name
    ])

    chromosome_order = sorted(combined_data['X.CHROM'].unique(), key=sort_comparison_chromosomes)
    

    figure = {
        'data': [
            go.Scatter(
                x=combined_data[combined_data['Mutation Type'] == 'Dense']['POS'],
                y=combined_data[combined_data['Mutation Type'] == 'Dense']['X.CHROM'],
                mode='markers',
                marker={
                    'size': 5,
                    'opacity': 0.8,
                    'symbol': 'circle',  
                    'color': 'grey', 
                },
                name='Dense',
            ),
            go.Scatter(
                x=combined_data[combined_data['Mutation Type'] == 'Sparse']['POS'],
                y=combined_data[combined_data['Mutation Type'] == 'Sparse']['X.CHROM'],
                mode='markers',
                marker={
                    'size': 5,
                    'opacity': 1.0,
                    'symbol': 'circle-open',  
                    'color': 'red', 
                },
                name='Sparse',
            ),
        ],
        'layout': {
            'title': title + ' - Mutations called by different surrogate parents',
            'xaxis': {'title': 'Chromosome Position'},
            'yaxis': {
                'title': 'Chromosome.SurrogateParent', 
                'categoryorder': 'array', 
                'categoryarray': chromosome_order,
                'tickfont': {'size': 5},
                'tickangle': 45 
            }
        }
    }

    return figure    
    
if __name__ == '__main__':
    app.run_server(host='0.0.0.0', port=8051, debug=True)

# Folder to save HTML files
output_folder = "compare_surrogate_tracts"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

for idx, title in enumerate(comparison_titles):
    fig = update_comparison_plot(idx)  
    plot(fig, filename=os.path.join(output_folder, f"{title.replace(' ', '_')}.html"), auto_open=False)

zip_file_name = "compare_surrogate_tracts_html_files.zip"

with zipfile.ZipFile(zip_file_name, 'w') as zipf:
    for root, _, files in os.walk(output_folder):
        for file in files:
            zipf.write(os.path.join(root, file), os.path.relpath(os.path.join(root, file), output_folder))

# Create a dictionary of data frames by name
data_frames_dict = {
    "C22asfather_C21askid": C22asfather_C21askid,
    "C23asfather_C21askid": C23asfather_C21askid,
    "C21asfather_C22askid": C21asfather_C22askid,
    "C23asfather_C22askid": C23asfather_C22askid,
    "C21asfather_C23askid": C21asfather_C23askid,
    "C22asfather_C23askid": C22asfather_C23askid,
    "C42asfather_C41askid": C42asfather_C41askid,
    "C41asfather_C42askid": C41asfather_C42askid,
    "P2P3asparents_P1askid": P2P3asparents_P1askid,
    "P2P4asparents_P1askid": P2P4asparents_P1askid,
    "P3P4asparents_P1askid": P3P4asparents_P1askid,
    "P1P3asparents_P2askid": P1P3asparents_P2askid,
    "P1P4asparents_P2askid": P1P4asparents_P2askid,
    "P3P4asparents_P2askid": P3P4asparents_P2askid,
    "P1P2asparents_P3askid": P1P2asparents_P3askid,
    "P1P4asparents_P3askid": P1P4asparents_P3askid,
    "P2P4asparents_P3askid": P2P4asparents_P3askid,
    "P1P2asparents_P4askid": P1P2asparents_P4askid,
    "P1P3asparents_P4askid": P1P3asparents_P4askid,
    "P2P3asparents_P4askid": P2P3asparents_P4askid,
}

comparison_groups = [
    ['C22asfather_C21askid', 'C23asfather_C21askid'],
    ['C21asfather_C22askid', 'C23asfather_C22askid'],
    ['C21asfather_C23askid', 'C22asfather_C23askid'],
    ['P2P3asparents_P1askid', 'P2P4asparents_P1askid', 'P3P4asparents_P1askid'],
    ['P1P3asparents_P2askid', 'P1P4asparents_P2askid', 'P3P4asparents_P2askid'],
    ['P1P4asparents_P3askid', 'P1P2asparents_P3askid', 'P2P4asparents_P3askid'],
    ['P1P2asparents_P4askid', 'P1P3asparents_P4askid', 'P2P3asparents_P4askid']
]

comparison_titles = [
    'C21_mutations',
    'C22_mutations',
    'C23_mutations',
    'P1_mutations',
    'P2_mutations',
    'P3_mutations',
    'P4_mutations'
]

def save_data_frames():
    # Create empty DataFrames for each comparison title
    data_frames = {
        title: pd.DataFrame(columns=['X.CHROM', 'POS', 'REF.ALT', 'Mutation Type', 'Source'])
        for title in comparison_titles
    }
    return data_frames

# Create a variable to store the data frames
mutations_data_frames = save_data_frames()

app = dash.Dash(__name__)

comparison_group_selector = dcc.Dropdown(
    id='comparison_group_selector',
    options=[{'label': title, 'value': i} for i, title in enumerate(comparison_titles)],
    value=0
)

app.layout = html.Div([
    comparison_group_selector,
    dcc.Graph(id='sparse-plot'),
    html.Div(id='output-container', style={"white-space": "pre-line"})
])

@app.callback(
    [Output('sparse-plot', 'figure'),
     Output('output-container', 'children')],
    [Input('comparison_group_selector', 'value')]
)


def update_plot(selected_index):
    global mutations_data_frames 
    selected_group = comparison_groups[selected_index]
    title = comparison_titles[selected_index]
    all_positions = {}
    traces = []
    mutation_counts = []

    all_mutations = pd.DataFrame(columns=['X.CHROM', 'POS', 'REF.ALT', 'Mutation Type', 'Source'])

    # Iterate through all the selected data frames and collect all the positions
    for frame_name in selected_group:
        mutation_data = data_frames_dict[frame_name]
        mutation_data = mutation_data[mutation_data['Mutation Type'] == 'Sparse']
        mutation_counts.append(f"{frame_name}: {len(mutation_data)}")


        for index, row in mutation_data.iterrows():
            pos = row['POS']
            chrom = row['X.CHROM']
            ref_alt_value = row['REF.ALT']
            if pos in all_positions:
                all_positions[pos]['count'] += 1
                all_positions[pos]['Sources'].append(frame_name)  
                all_positions[pos]['REF.ALT'].append(ref_alt_value)
            else:
                all_positions[pos] = {'count': 1, 'chrom': chrom, 'Sources': [frame_name], 'REF.ALT': [ref_alt_value]}  

    x_values = []
    y_values = []
    colors = []
    hover_texts = []  # Use hover_texts instead of frame_names

    # Iterate through the collected positions and create the trace
    for pos, info in all_positions.items():
        if info['count'] > 1:
            color = 'purple' # Color for overlapping mutations
            hover_text = 'Shared' # If shared
        else:
            color = 'orange' # Color for non-overlapping mutations
            hover_text = info['Sources'][0] # Name of the unique data frame

        x_values.append(pos)
        y_values.append(info['chrom'])
        colors.append(color)
        hover_texts.append(hover_text) # Append hover text
            
    # Iterate through the collected positions and update DataFrames
    for pos, info in all_positions.items():
        source_value = "Shared" if info['count'] > 1 else info['Sources'][0]
        ref_alt_value = ','.join(set(info['REF.ALT']))
        new_row = pd.DataFrame({
            'X.CHROM': [info['chrom']],
            'POS': [pos],
            'REF.ALT': [ref_alt_value], 
            'Mutation Type': ['Sparse'],
            'Source': [source_value]
        })
        all_mutations = pd.concat([all_mutations, new_row], ignore_index=True)

    # Update the corresponding DataFrame in mutations_data_frames
    mutations_data_frames[title] = all_mutations
    
    chromosome_order = sorted(list(set(y_values)), key=sort_chromosomes)

    trace = go.Scatter(
        x=x_values,
        y=y_values,
        mode='markers',
        marker={
            'size': 10,
            'opacity': 1.0,
            'symbol': 'circle-open',
            'color': colors
        },
        text=hover_texts, # Use hover_texts for text attribute
    )

    traces.append(trace)

    figure = {
        'data': traces,
        'layout': {
            'title': title,
            'xaxis': {'title': 'Chromosome Position'},
            'yaxis': {'title': 'Chromosome', 'categoryorder': 'array', 'categoryarray': chromosome_order}
        }
    }
    
    unique_count = len(all_mutations[all_mutations['Source'] != "Shared"])
    shared_count = len(all_mutations[all_mutations['Source'] == "Shared"])

    output_text = f"\n\nUnique Mutations: {unique_count}"
    output_text += f"\nShared Mutations: {shared_count}"
    
    return figure, output_text

if __name__ == '__main__':
    app.run_server(host='0.0.0.0', port=8052, debug=True)

# Folder to save HTML files
output_folder = "shared_and_unique_html_files"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

for idx, title in enumerate(comparison_titles):
    fig, _ = update_plot(idx)  # You may need to adjust this part based on how you generate your figures
    plot(fig, filename=os.path.join(output_folder, f"{title.replace(' ', '_')}.html"), auto_open=False)

zip_file_name = "shared_and_unique_html_files.zip"

with zipfile.ZipFile(zip_file_name, 'w') as zipf:
    for root, _, files in os.walk(output_folder):
        for file in files:
            zipf.write(os.path.join(root, file), os.path.relpath(os.path.join(root, file), output_folder))

# Filter for sparse mutations
C41_mutations_filtered = data_frames_dict["C42asfather_C41askid"]
C41_mutations_filtered = C41_mutations_filtered[C41_mutations_filtered['Mutation Type'] == 'Sparse']

C42_mutations_filtered = data_frames_dict["C41asfather_C42askid"]
C42_mutations_filtered = C42_mutations_filtered[C42_mutations_filtered['Mutation Type'] == 'Sparse']


# Define sparse shared/unique mutations as data frames for downstream analysis
C21_mutations = mutations_data_frames['C21_mutations']
C22_mutations = mutations_data_frames['C22_mutations']
C23_mutations = mutations_data_frames['C23_mutations']
P1_mutations = mutations_data_frames['P1_mutations']
P2_mutations = mutations_data_frames['P2_mutations']
P3_mutations = mutations_data_frames['P3_mutations']
P4_mutations = mutations_data_frames['P4_mutations']

app = dash.Dash(__name__)

# Define the layout
app.layout = html.Div([
    dcc.Dropdown(
        id='data-selector',
        options=[
            {'label': 'C11 Mutations', 'value': 'C11'},
            {'label': 'C12 Mutations', 'value': 'C12'},
            {'label': 'C21 Mutations', 'value': 'C21'},
            {'label': 'C22 Mutations', 'value': 'C22'},
            {'label': 'C23 Mutations', 'value': 'C23'},
            {'label': 'C31 Mutations', 'value': 'C31'},
            {'label': 'C32 Mutations', 'value': 'C32'},
            {'label': 'C41 Mutations', 'value': 'C41'},
            {'label': 'C42 Mutations', 'value': 'C42'},
            {'label': 'P1 Mutations', 'value': 'P1'},
            {'label': 'P2 Mutations', 'value': 'P2'},
            {'label': 'P3 Mutations', 'value': 'P3'},
            {'label': 'P4 Mutations', 'value': 'P4'}
        ],
        value='C21' # Default selection
    ),
    dcc.Graph(id='mutation-plot')
])

# Define callback to update plot
@app.callback(
    Output('mutation-plot', 'figure'),
    [Input('data-selector', 'value')]
)
def update_plot(selected_data):
    # Map selected value to corresponding DataFrame
    data_mapping = {
        'C11': C11askid,
        'C12': C12askid,
        'C21': C21_mutations,
        'C22': C22_mutations,
        'C23': C23_mutations,
        'C31': C31askid,
        'C32': C32askid,
        'C41': C41_mutations_filtered,
        'C42': C42_mutations_filtered,
        'P1': P1_mutations,
        'P2': P2_mutations,
        'P3': P3_mutations,
        'P4': P4_mutations,
    }
    
    # Create the scatter plot for the selected DataFrame
    fig = px.scatter(data_mapping[selected_data], x='POS', y='X.CHROM', title=f'{selected_data} Mutations')

    return fig

# Run the app
if __name__ == '__main__':
    app.run_server(host='0.0.0.0', port=8053, debug=True)



directory_path = "final_denoms"

# Initialize the total accessible bases count
total_accessible_bases_count = 0

for filename in os.listdir(directory_path):
    if filename.endswith(".bed"):
        file_path = os.path.join(directory_path, filename)

        accessible_df = pd.read_csv(file_path, sep="\t", header=None, names=["Chromosome", "Start", "End"])

        accessible_bases = accessible_df["End"].sum() - accessible_df["Start"].sum()
        
        total_accessible_bases_count += accessible_bases

        print(f"Processed file {filename}. Accessible bases: {accessible_bases}")

data_frames = {
    'C11': C11askid,
    'C12': C12askid,
    'C21': C21_mutations,
    'C22': C22_mutations,
    'C23': C23_mutations,
    'C31': C31askid,
    'C32': C32askid,
    'C41': C41_mutations_filtered,
    'C42': C42_mutations_filtered,
    'P1': P1_mutations,
    'P2': P2_mutations,
    'P3': P3_mutations,
    'P4': P4_mutations,
}


denominators = {
    'C11': 2701991570,
    'C12': 2701991570,
    'C21': 1919372907,
    'C22': 2091899533,
    'C23': 2141739396,
    'C31': 2701991570,
    'C32': 2701991570,
    'C41': 1164137259,
    'C42': 1169041661,
    'P1': 1868756884,
    'P2': 1810864496,
    'P3': 1471216193,
    'P4': 1779830179,
}

# Create a list of data frame names and their corresponding counts
data_frame_names = list(data_frames.keys())

# Bar plot 1: Mutation counts
mutation_counts = [len(df) for df in data_frames.values()]

# Bar plot 2: Denominators * 2
denominator_times_2 = [denominators[name] * 2 for name in data_frame_names]

# Bar plot 3: Mutation counts divided by denominators * 2
normalized_counts = [count / (denominators[name] * 2) for name, count in zip(data_frame_names, mutation_counts)]

fig, axs = plt.subplots(3, 1, figsize=(10, 15))

axs[0].bar(data_frame_names, mutation_counts, color='b', alpha=0.7)
axs[1].bar(data_frame_names, denominator_times_2, color='g', alpha=0.7)
axs[2].bar(data_frame_names, normalized_counts, color='r', alpha=0.7)

axs[0].set_title('Mutation Counts by Data Frame')
axs[0].set_ylabel('Mutation Count')
axs[1].set_title('Denominators * 2 by Data Frame')
axs[1].set_ylabel('Variance in Accessible Bases*2')
axs[2].set_title('Normalized Mutation Counts by Data Frame')
axs[2].set_ylabel('Mutation Rates')

formatter = ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((-1,1))
axs[1].yaxis.set_major_formatter(formatter)
axs[2].yaxis.set_major_formatter(formatter)

axs[2].axhline(y=1.29*10**-8, color='black', linestyle='--', linewidth=1)

# Setting the x-axis label for all subplots
for ax in axs:
    ax.set_xlabel('Individual')

plt.tight_layout()
plt.show()

tally_list = []

mut_types = ['A.C', 'A.G', 'A.T', 'C.T', 'C.G', 'C.A']

for child_id, df in data_frames.items():
    if 'REF.ALT' in df.columns:
        tally = df['REF.ALT'].value_counts().reset_index()
        tally.columns = ['mut_type', 'count']
        all_mut_types = pd.DataFrame({'mut_type': mut_types}) 
        tally = pd.merge(all_mut_types, tally, on='mut_type', how='left').fillna(0)
        tally['child_ID'] = child_id
        tally = tally[['child_ID', 'mut_type', 'count']]
        tally_list.append(tally)

# Concatenate all the tally DataFrames
final_tally = pd.concat(tally_list, ignore_index=True)

final_tally.to_csv("final_tally.csv", index=False)

for name, df in data_frames.items():
        # Create the BED formatted dataframe
    bed_df = pd.DataFrame({
        'chrom': df['X.CHROM'],
        'chromStart': df['POS'] - 1,
        'chromEnd': df['POS']
    })
    
    # Save the dataframe to a BED file
    bed_df.to_csv(f'{name}.bed', sep='\t', index=False, header=False)
