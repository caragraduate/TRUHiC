import pandas as pd
import os

def expand_coordinates(row):
    """Expand loop coordinates by ±5kb."""
    return pd.Series({
        'min_x': row['x1'] - 5000,
        'max_x': row['x2'] + 5000,
        'min_y': row['y1'] - 5000,
        'max_y': row['y2'] + 5000
    })

def calculate_center_and_area(row):
    """Calculate the center and expand the search area by ±5kb."""
    center_x = (row['x1'] + row['x2']) / 2
    center_y = (row['y1'] + row['y2']) / 2
    return pd.Series({
        'center_x': center_x, 'center_y': center_y,
        'min_x': center_x - 5000, 'max_x': center_x + 5000,
        'min_y': center_y - 5000, 'max_y': center_y + 5000
    })

def find_overlapping_loops(row1, df2):
    """Find overlapping loops within tolerance range."""
    overlaps = df2[
        (df2['#chr1'] == row1['#chr1']) &
        (df2['min_x'] <= row1['max_x']) & (df2['max_x'] >= row1['min_x']) &
        (df2['min_y'] <= row1['max_y']) & (df2['max_y'] >= row1['min_y'])
    ]
    return overlaps.index.tolist() if not overlaps.empty else []

def calculate_jaccard_index(enhanced_loop_callset, HR_loop_callset):
    """Calculate the Jaccard Index between enhanced loops and HR loops."""
    # Find overlaps for each loop
    enhanced_loop_callset['overlapping_loops_in_HR'] = enhanced_loop_callset.apply(
        lambda row: find_overlapping_loops(row, HR_loop_callset), axis=1
    )
    HR_loop_callset['overlapping_loops_in_enhanced'] = HR_loop_callset.apply(
        lambda row: find_overlapping_loops(row, enhanced_loop_callset), axis=1
    )

    # Count the number of loops with at least one overlap
    overlaps_in_enhanced = (enhanced_loop_callset['overlapping_loops_in_HR'].apply(len) > 0).sum()
    overlaps_in_HR = (HR_loop_callset['overlapping_loops_in_enhanced'].apply(len) > 0).sum()
    # Print results
    print(f"Number of overlapping loops in loop1: {overlaps_in_enhanced}")
    print(f"Number of overlapping loops in loop2: {overlaps_in_HR}")

    # Collect all unique overlapping loop indices
    def find_total_overlaps(df):
        all_overlaps = set()
        for overlaps in df['overlapping_loops_in_HR']:
            all_overlaps.update(overlaps)  # Add all unique overlapping indices
        return len(all_overlaps)

    # Calculate the total unique overlapping loops
    total_unique_overlaps = find_total_overlaps(enhanced_loop_callset)

    # Calculate the union of loops
    total_loops_enhanced = len(enhanced_loop_callset)
    total_loops_HR = len(HR_loop_callset)
    union_of_loops = total_loops_enhanced + total_loops_HR - total_unique_overlaps

    # Calculate Jaccard Index
    jaccard_index = total_unique_overlaps / union_of_loops if union_of_loops > 0 else 0

    # Print the Jaccard Index
    print("Number of Overlapping Loops:", total_unique_overlaps)
    print("Number of Total Loops:", union_of_loops)
    print(f"Jaccard Index: {jaccard_index}")

    return round(jaccard_index, 4)

def calculate_metrics(loop_file, HR_loop_file, chrom):
    """Calculate metrics (F1 score, overlaps, Jaccard index) for a given chromosome."""
    try:
        # Load enhanced loop file
        enhanced_loop = os.path.join(loop_file, f"hiccups_results_{chrom}/merged_loops.bedpe")
        enhanced_loop_callset = pd.read_csv(enhanced_loop, sep='\t', skiprows=[1])
        # Sorting by chromosome and start coordinate
        enhanced_loop_callset.sort_values(by=['#chr1', 'x1', 'y1'], inplace=True)
        enhanced_loop_callset.reset_index(drop=True, inplace=True)

        # Load HR loop file
        # for human
        HR_loop = os.path.join(HR_loop_file, f"hiccups_results_ori_KR_{chrom}/merged_loops.bedpe")
        # # for mouse: 
        # HR_loop = os.path.join(HR_loop_file, f"hiccups_results_ori_KR_{chrom}/merged_loops.bedpe")
        HR_loop_callset = pd.read_csv(HR_loop, sep='\t', skiprows=[1])
        # Sorting by chromosome and start coordinate
        HR_loop_callset.sort_values(by=['#chr1', 'x1', 'y1'], inplace=True)
        HR_loop_callset.reset_index(drop=True, inplace=True)

        if enhanced_loop_callset.empty or HR_loop_callset.empty:
            print(f"Skipping {chrom} due to empty data in files.")
            return None
    except FileNotFoundError as e:
        print(f"File not found for {chrom}: {e}")
        return None

    # Expand coordinates and calculate center
    enhanced_loop_callset[['min_x', 'max_x', 'min_y', 'max_y']] = enhanced_loop_callset.apply(expand_coordinates, axis=1)
    HR_loop_callset[['min_x', 'max_x', 'min_y', 'max_y']] = HR_loop_callset.apply(expand_coordinates, axis=1)

    # Count loops
    enhanced_loop_count = len(enhanced_loop_callset)
    HR_loop_count = len(HR_loop_callset)

    # Calculate overlaps
    enhanced_loop_callset['overlapping_loops_in_HR'] = enhanced_loop_callset.apply(
        lambda row: find_overlapping_loops(row, HR_loop_callset), axis=1
    )
    HR_loop_callset['overlapping_loops_in_enhanced'] = HR_loop_callset.apply(
        lambda row: find_overlapping_loops(row, enhanced_loop_callset), axis=1
    )

    TP = (enhanced_loop_callset['overlapping_loops_in_HR'].apply(len) > 0).sum()
    FP = (enhanced_loop_callset['overlapping_loops_in_HR'].apply(len) == 0).sum()
    FN = (HR_loop_callset['overlapping_loops_in_enhanced'].apply(len) == 0).sum()

    # Calculate F1 score
    F1_score = TP / (TP + 0.5 * (FP + FN)) if (TP + 0.5 * (FP + FN)) > 0 else 0
    
    # Apply the function to both DataFrames
    enhanced_loop_callset[['center_x', 'center_y', 'min_x', 'max_x', 'min_y', 'max_y']] = enhanced_loop_callset.apply(calculate_center_and_area, axis=1)
    HR_loop_callset[['center_x', 'center_y', 'min_x', 'max_x', 'min_y', 'max_y']] = HR_loop_callset.apply(calculate_center_and_area, axis=1)

    # Calculate Jaccard Index
    jaccard_index = calculate_jaccard_index(enhanced_loop_callset, HR_loop_callset)

    return TP, FP, FN, F1_score, jaccard_index, enhanced_loop_count, HR_loop_count

def main():
    base_dir = "/your/path/to/base_dir"
    replicates = ["rep4"]
    base_model_dirs = [
        "/your/path/to/results",
        # Add more folders here as needed
    ]
    chromosomes = ['chr18', 'chr19', 'chr20', 'chr21', 'chr22']
    HR_loop_file = '/Results/GM12878'
    
    output_results = []

    # Generate directories dynamically for all replicates
    model_dirs = []
    for replicate in replicates:
        for base_model_dir in base_model_dirs:
            # updated_dir = base_model_dir.replace("rep1", "rep1").replace("rep_1", f"{replicate.replace('rep', 'rep_')}")
            updated_dir = base_model_dir.replace("rep1", replicate)
            model_dirs.append(updated_dir)

    for model_dir in model_dirs:
        replicate = model_dir.split("_")[-2]  # Extract replicate dynamically (e.g., rep1, rep2, rep3)
        model_name = model_dir  # Extract model name dynamically

        print(f"Processing Model: {model_name}, Replicate: {replicate}")
        loop_file = os.path.join(base_dir, model_dir)

        for chrom in chromosomes:
            print(f"  Processing Chromosome: {chrom}")
            result = calculate_metrics(loop_file, HR_loop_file, chrom)
            if result is None:
                continue

            TP, FP, FN, F1_score, jaccard_index, enhanced_loop_count, HR_loop_count = result
            output_results.append({
                "Model": model_name,
                "Replicate": replicate,
                "Chromosome": chrom,
                "HR Loop Count": HR_loop_count,
                "True Positives (TP)": TP,
                "False Positives (FP)": FP,
                "False Negatives (FN)": FN,
                "Enhanced Loop Count": enhanced_loop_count,
                "Jaccard Index": jaccard_index,
                "F1 Score": round(F1_score, 4),

            })

    # Save results to CSV
    output_df = pd.DataFrame(output_results)
    output_csv = "Loop_F1_scores_Jaccard_results_benchmark.csv"
    output_df.to_csv(output_csv, index=False, sep='\t')
    print(f"Results saved to {output_csv}")

if __name__ == "__main__":
    main()