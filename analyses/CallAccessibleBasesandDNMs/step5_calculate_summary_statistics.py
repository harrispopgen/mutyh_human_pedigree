import os
import pandas as pd

# Read the average depth files and calculate summary statistics
output_directory = "/net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/samtools_depth_files/averaged_depths_per_chr/output"

summary_stats = []

# Initialize the variable for counting accessible bases
accessible_bases_count = 0

for filename in os.listdir(output_directory):
    if filename.endswith(".txt"):
        file_path = os.path.join(output_directory, filename)
        chromosome_df = pd.read_csv(file_path, sep="\t", header=None, names=["Chromosome", "Start", "End", "Avg_Depth"])

        chromosome_name = chromosome_df['Chromosome'].iloc[0]
        mean_depth = chromosome_df['Avg_Depth'].mean()
        median_depth = chromosome_df['Avg_Depth'].median()
        min_depth = chromosome_df['Avg_Depth'].min()
        max_depth = chromosome_df['Avg_Depth'].max()
        mean_2_5 = 2.5 * mean_depth

        # Calculate the number of accessible bases
        accessible_bases_file = "/net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/samtools_depth_files/averaged_depths_per_chr/output/averaged_accessible_bases.bed"
        accessible_df = pd.read_csv(accessible_bases_file, sep=",")  
        accessible_df.columns = ["Chromosome", "Start", "End", "Avg_Depth"]  
        accessible_bases = accessible_df.loc[accessible_df["Chromosome"] == chromosome_name, "End"].sum() - accessible_df.loc[accessible_df["Chromosome"] == chromosome_name, "Start"].sum()
        accessible_bases_count += accessible_bases

        summary_stats.append({
            'Chromosome': chromosome_name,
            'Mean': mean_depth,
            'Median': median_depth,
            'Range': f"{min_depth} - {max_depth}",
            '2.5xMean': mean_2_5
        })

# Create summary table DataFrame
summary_table = pd.DataFrame(summary_stats)

# Print the summary table
print(summary_table)

# Save the summary table as a CSV file
summary_table.to_csv("/net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/samtools_depth_files/averaged_depths_per_chr/output/summary_table.csv", index=False)

# Print the number of accessible bases
print(f"Number of accessible bases: {accessible_bases_count}")

# Convert averaged_accessible_bases.bed to tab-separated .bed file
accessible_bed_file = "/net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/samtools_depth_files/averaged_depths_per_chr/output/averaged_accessible_bases.bed"
output_bed_file = "/net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/samtools_depth_files/averaged_depths_per_chr/output/averaged_accessible_bases_no_depths.bed"

# Read the bed file and remove the depth column
accessible_df = pd.read_csv(accessible_bed_file, sep=",")
accessible_df.drop(columns="Avg_Depth", inplace=True)

# Save the modified bed file as tab-separated .bed
accessible_df.to_csv(output_bed_file, sep="\t", header=False, index=False)
