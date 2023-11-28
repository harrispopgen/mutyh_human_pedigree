import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import matplotlib.patches as patches

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

