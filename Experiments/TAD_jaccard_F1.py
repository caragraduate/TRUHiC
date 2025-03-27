import pandas as pd
import os
import subprocess

def process_tads(TAD_file, chrom, ratio):
    """Process TADs for a given TAD file and chromosome."""
    # Load the enhanced TAD file
    enhanced_TAD = os.path.join(TAD_file, f'preds_lr_test_{chrom}_ratio{ratio}_convert_10kb/10000_blocks.bedpe')
    TAD_10kb = pd.read_csv(enhanced_TAD, sep='\t', skiprows=[1])

    # Extract relevant columns and sort
    TAD_10kbnew = TAD_10kb[['#chr1', 'x1', 'x2']]
    TAD_10kbnew.columns = ['chr', 'TAD_start', 'TAD_end']
    TAD_10kbnew_1 = TAD_10kbnew.sort_values(['chr', 'TAD_start'], ignore_index=True)

    # Save sorted TADs
    output_file = os.path.join(TAD_file, f'preds_lr_test_{chrom}_ratio{ratio}_convert_10kb/{chrom}_TADs_ratio{ratio}.bed')
    TAD_10kbnew_1.to_csv(output_file, index=False, sep='\t', header=True)

    return len(TAD_10kbnew_1), output_file


def calculate_jaccard(file_a, file_b):
    """Calculate Jaccard index using bedtools."""
    # Sort only the first file (enhanced TAD file)
    sorted_file_a = file_a.replace(".bed", "_sorted.bed")
    subprocess.run(f"sort -k1,1 -k2,2n {file_a} > {sorted_file_a}", shell=True, check=True)

    # Calculate Jaccard index using bedtools
    jaccard_result = subprocess.run(
        f"bedtools jaccard -a {sorted_file_a} -b {file_b}",
        shell=True, capture_output=True, text=True
    )
    if jaccard_result.returncode == 0:
        # Extract the Jaccard index (third column from the output)
        lines = jaccard_result.stdout.split("\n")
        if len(lines) > 1:  # Ensure the result has at least two lines (header + values)
            values = lines[1].split("\t")  # Split the second line (results) by tab
            jaccard_index = float(values[2])  # Third column is the Jaccard value
            return jaccard_index
    else:
        print(f"Error calculating Jaccard index: {jaccard_result.stderr}")
        return None


def calculate_f1_score(enhanced_file, HR_file):
    """Calculate F1 score for a given enhanced and HR TAD file."""
    # Load TAD files
    enhanced_tad_callset = pd.read_csv(enhanced_file, sep='\t', header=None, skiprows=1)
    HR_tad_callset = pd.read_csv(HR_file, sep='\t', header=None, skiprows=1)

    # Function to check overlaps between two TAD regions
    def find_overlapping_tads(row1, df2):
        # Ensure numeric type conversion locally for relevant columns
        df2[1] = df2[1].astype(int)
        df2[2] = df2[2].astype(int)
        row1[1] = int(row1[1])
        row1[2] = int(row1[2])

        # Find overlapping TADs
        overlaps = df2[(df2[0] == row1[0]) &  # Compare chromosome
                    (df2[1] <= row1[2]) &  # Check if start in df2 <= end in row1
                    (df2[2] >= row1[1])]   # Check if end in df2 >= start in row1
        return overlaps.index.tolist() if not overlaps.empty else []

    # Find overlaps for each dataset
    enhanced_tad_callset['overlapping_tads_in_HR'] = enhanced_tad_callset.apply(
        lambda row: find_overlapping_tads(row, HR_tad_callset), axis=1
    )
    HR_tad_callset['overlapping_tads_in_enhanced'] = HR_tad_callset.apply(
        lambda row: find_overlapping_tads(row, enhanced_tad_callset), axis=1
    )

    # Calculate metrics
    TP = (enhanced_tad_callset['overlapping_tads_in_HR'].apply(len) > 0).sum()  # True Positives
    FP = (enhanced_tad_callset['overlapping_tads_in_HR'].apply(len) == 0).sum()  # False Positives
    FN = (HR_tad_callset['overlapping_tads_in_enhanced'].apply(len) == 0).sum()  # False Negatives

    # Calculate F1 score
    if (TP + 0.5 * (FP + FN)) > 0:
        F1_score = TP / (TP + 0.5 * (FP + FN))
    else:
        F1_score = 0

    return round(F1_score, 4)

def main():
    base_dir = "/your/path/to/base_dir"
    # replicates = ["rep1", "rep2", "rep3"]
    replicates = ["rep4"]
    base_model_dirs = [
        "/your/path/to/results",
        # Add more folders here as needed
    ]
    chromosomes = ['chr18', 'chr19', 'chr20', 'chr21', 'chr22']
    ratio = 16
    HR_base_dir = '/Results/GM12878'

    output_results = []

    # Generate directories dynamically for all replicates
    model_dirs = []
    for replicate in replicates:
        for base_model_dir in base_model_dirs:
            # Replace "rep1" in the directory name with the current replicate
            updated_dir = base_model_dir.replace("rep1", "rep1")
            # Append the directory to the base path
            model_dirs.append(os.path.join(base_dir, updated_dir))

    for model_dir in model_dirs:
        # Extract replicate dynamically from the directory name
        replicate = model_dir.split("_")[-2]  # Adjust index if needed
        # Extract model name from the directory path
        model_name = model_dir  # Extract model name dynamically
        TAD_file = model_dir

        print(f"Processing Model: {model_name}, Replicate: {replicate}")

        for chrom in chromosomes:
            print(f"  Processing Chromosome: {chrom}")
            
            # Process TAD boundaries and get the detected number
            detected_count, enhanced_TAD_file = process_tads(TAD_file, chrom, ratio)

            # Define the HR TAD boundary file
            HR_TAD_file = os.path.join(HR_base_dir, f'GM12878_ori_{chrom}_10kb/HR_{chrom}_TADs_ratio16.bedpe')

            # Calculate Jaccard index
            jaccard_index = calculate_jaccard(enhanced_TAD_file, HR_TAD_file)

            # Calculate F1 score
            F1_score = calculate_f1_score(enhanced_TAD_file, HR_TAD_file)

            # Append results
            output_results.append({
                "Model": model_name,
                "Replicate": replicate,
                "Chromosome": chrom,
                "Detected TADs": detected_count,
                "Jaccard Index": round(jaccard_index, 4) if jaccard_index is not None else None,
                "F1 Score": F1_score
            })

    # Save results to CSV
    output_df = pd.DataFrame(output_results)
    output_csv = "TAD_Jaccard_F1_results_benchmark_GM12878_Feb05.csv"
    output_df.to_csv(output_csv, index=False, sep='\t')
    print(f"Results saved to {output_csv}")


if __name__ == "__main__":
    main()