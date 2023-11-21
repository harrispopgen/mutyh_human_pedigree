import os
import pandas as pd
import random
import numpy as np

# Set sampling percentage
sampling_percentage = 0.001

# Initialize with random seed for reproducibility
random.seed(42)

# Sample IDs and chromosomes
sample_ids = ["C11.1079256", "C12.1079257", "C21.1079259", "C22.1079260", "C23.1079261", "C31.1079264", "C32.1079265", "C41.1079267", "C42.1079268", "P1.1079254", "P2.1079258", "P3.1079262", "P4.1079266", "S1.1079255", "S3.1079263"]
chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"]

error_directory = "/net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/samtools_depth_files/averaged_depths_per_chr/errors"
if not os.path.exists(error_directory):
    os.makedirs(error_directory)

def log_error(message):
    with open(os.path.join(error_directory, "error.log"), "a") as error_file:
        error_file.write(message + "\n")

# Function to merge data for all samples and chromosomes
def merge_data():
    try:
        # Check that txt files per chromosome across individuals have the same length
        chromosome_counts = pd.DataFrame(columns=["sample_id", "chromosome", "num_windows"])

        for sample_id in sample_ids:
            for chromosome in chromosomes:
                file_name = f"{sample_id}_{chromosome}.txt"
                file_path = os.path.join("/net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/samtools_depth_files/depths_per_chr_and_individual/output", file_name)
                if os.path.exists(file_path):
                    data = pd.read_csv(file_path, sep="\t", header=None, names=["Chromosome", "Start", "End", "Depth"])
                    num_windows = len(data)
                    chromosome_counts = chromosome_counts.append({"sample_id": sample_id, "chromosome": chromosome, "num_windows": num_windows}, ignore_index=True)

        # Group the data frame by chromosome
        grouped_counts = chromosome_counts.groupby("chromosome")["num_windows"].nunique().reset_index()
        grouped_counts.rename(columns={"num_windows": "same_num_windows"}, inplace=True)

        print(grouped_counts)

        # Create merged data frame of all depth scores per 10KB window, per chromosome across individuals
        output_data = pd.DataFrame(columns=["Chromosome", "Start", "End", "Avg_Depth"] + sample_ids)

        for chromosome in chromosomes:
            # Read the input file for the first sample to determine the number of rows
            sample_id = sample_ids[0]
            file_name = f"{sample_id}_{chromosome}.txt"
            file_path = os.path.join("/net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/samtools_depth_files/depths_per_chr_and_individual/output", file_name)

            if os.path.exists(file_path):
                sample_file_data = pd.read_csv(file_path, sep="\t", header=None, names=["Chromosome", "Start", "End", "Depth"])
                n_rows = len(sample_file_data)
            else:
                print(f"File {file_name} does not exist.")
                continue

            # Initialize an empty data frame to store sample data for the current chromosome
            sample_data = pd.DataFrame(index=range(n_rows), columns=sample_ids)

            for sample_id in sample_ids:
                file_name = f"{sample_id}_{chromosome}.txt"
                file_path = os.path.join("/net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/samtools_depth_files/depths_per_chr_and_individual/output", file_name)

                if os.path.exists(file_path):
                    sample_file_data = pd.read_csv(file_path, sep="\t", header=None, names=["Chromosome", "Start", "End", sample_id])
                    sample_data[sample_id] = sample_file_data[sample_id]
                else:
                    print(f"File {file_name} does not exist.")

            # Add the chromosome information to the sample_data data frame
            sample_data["Chromosome"] = chromosome
            sample_data["Start"] = sample_file_data["Start"]
            sample_data["End"] = sample_file_data["End"]

            # Append the chromosome data to the output_data data frame
            output_data = pd.concat([output_data, sample_data], ignore_index=True)

        # Check for any duplicates in the merged data frame
        duplicates = output_data[output_data.duplicated(subset=["Chromosome", "Start", "End"], keep=False)]
        print("Potential Duplicates:")
        print(duplicates)

        # Check the number of rows for each chromosome and sample
        print("Number of rows for each sample and chromosome:")
        grouped_output_data = output_data.groupby(["Chromosome"] + sample_ids).size().reset_index(name="NumRows")
        print(grouped_output_data)

        print(output_data)

        # Calculate the average depth score per row
        output_data['Avg_Depth'] = output_data[sample_ids].mean(axis=1)

        # Create separate data frames for each chromosome containing required columns
        chromosome_data = output_data[['Chromosome', 'Start', 'End', 'Avg_Depth']].groupby('Chromosome')

        # Save separate data frames as text files
        output_directory = "/net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/samtools_depth_files/averaged_depths_per_chr/output"
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)

        for chromosome, df in chromosome_data:
            chromosome_name = df['Chromosome'].iloc[0]
            file_name = f"{chromosome_name}_average_depth.txt"
            file_path = os.path.join(output_directory, file_name)
            df.to_csv(file_path, sep="\t", index=False, header=False)

        # Check output data matches data in original files
        # Calculate the number of sites to sample
        total_sites = len(pd.unique(output_data["Chromosome"])) * len(pd.unique(output_data["Start"]))
        sites_to_sample = round(total_sites * (sampling_percentage / 100))

        # Randomly select sites to sample
        sampled_sites = random.sample(range(1, total_sites + 1), sites_to_sample)

        # Initialize counters
        total_comparisons = 0
        matched_cases = 0

        # Iterate over the sampled sites
        for site_index in sampled_sites:
            # Get the chromosome and start position for the current site
            chromosome = pd.unique(output_data["Chromosome"])[(site_index - 1) // len(pd.unique(output_data["Start"]))]
            start = pd.unique(output_data["Start"])[(site_index - 1) % len(pd.unique(output_data["Start"]))]

            # Iterate over the sample IDs
            for sample_id in sample_ids:
                file_name = f"{sample_id}_{chromosome}.txt"
                file_path = os.path.join("/net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/samtools_depth_files/depths_per_chr_and_individual/output", file_name)

                if os.path.exists(file_path):
                    original_data = pd.read_csv(file_path, sep="\t", header=None, names=["Chromosome", "Start", "End", sample_id])

                    # Check if there are matching rows in the original_data DataFrame
                    matching_rows = original_data.loc[
                        (original_data["Chromosome"] == chromosome) & (original_data["Start"] == start), sample_id
                    ]

                    if not matching_rows.empty:
                        # Retrieve the value from the original text file
                        original_value = matching_rows.values[0]

                        # Retrieve the value from output_data
                        output_value = output_data.loc[
                            (output_data["Chromosome"] == chromosome) & (output_data["Start"] == start), sample_id
                        ].values[0]

                        total_comparisons += 1

                        if original_value == output_value:
                            matched_cases += 1
                            print(f"The values for site {chromosome} {start} and sample {sample_id} match.")
                        else:
                            print(f"The values for site {chromosome} {start} and sample {sample_id} do not match.")
                    else:
                        continue
                else:
                    print(f"File {file_name} does not exist.")

        percentage_matched = (matched_cases / total_comparisons) * 100
        print(f"Percentage of matched cases: {percentage_matched:.2f}%")

    except Exception as e:
        log_error(f"An error occurred during the merge_data function: {str(e)}")

# Calculate average depth per row and filter data frame
def calculate_average_depth_and_filter():
    try:
        # Read average depth files
        input_directory = "/net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/samtools_depth_files/averaged_depths_per_chr/output"
        averaged_depth_files = os.listdir(input_directory)

        # Initialize empty data frame
        filtered_data = pd.DataFrame(columns=["Chromosome", "Start", "End", "Avg_Depth"])

        # Iterate over files and filter data
        for file_name in averaged_depth_files:
            file_path = os.path.join(input_directory, file_name)
            if os.path.isfile(file_path):
                data = pd.read_csv(file_path, sep="\t", header=None, names=["Chromosome", "Start", "End", "Avg_Depth"])

                # Apply depth filters
                filtered_data = pd.concat([filtered_data, data[(data["Avg_Depth"] > 12) & (data["Avg_Depth"] < 120)]], ignore_index=True)

        # Save filtered data frame to new file
        output_file = "/net/harris/vol1/home/clyoung1/myh_pedigree/230705_fullpipeline/samtools_depth_files/averaged_depths_per_chr/output/averaged_accessible_bases.bed"
        filtered_data.to_csv(output_file, index=False)

    except Exception as e:
        log_error(f"An error occurred during the calculate_average_depth_and_filter function: {str(e)}")

# Execute functions
merge_data()
calculate_average_depth_and_filter()

