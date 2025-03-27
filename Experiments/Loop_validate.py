import pandas as pd
import os

# ==============================
# User-defined Parameters
# ==============================

cell_line = "GM12878"  # Change this for each cell line
models = [
    "/your/path/to/results",
    # Add more folders here as needed
]
chromosomes = ["chr18", "chr19", "chr20", "chr21", "chr22"]

# ==============================
# Helper Functions
# ==============================

def convert_chr_to_int(df, filename="Unknown File"):
    """ Convert chromosome column to integers, removing 'chr' prefix. """
    df['chr'] = df['chr'].astype(str).str.replace('chr', '', regex=True)
    df['chr'] = pd.to_numeric(df['chr'], errors='coerce')

    missing_chr_rows = df[df['chr'].isna()]
    if not missing_chr_rows.empty:
        print(f"⚠ Warning: Dropping {len(missing_chr_rows)} rows with non-numeric chromosomes in {filename}!")

    return df.dropna().astype({'chr': 'int'})

def expand_and_merge_intervals(df):
    """ Expands intervals by ±5kb and merges overlapping or adjacent intervals. """
    if df.empty:
        print("⚠ Warning: No loop intervals found to expand!")
        return pd.DataFrame(columns=['chr', 'start', 'end'])

    df['start'] -= 5000  # Expand start by -5kb
    df['end'] += 5000  # Expand end by +5kb
    df = df.sort_values(by=['chr', 'start'])

    merged_intervals = []
    current_chr, current_start, current_end = df.iloc[0]

    for i in range(1, len(df)):
        row_chr, row_start, row_end = df.iloc[i]

        if row_chr == current_chr and row_start <= current_end + 1:  # Merge overlapping/adjacent
            current_end = max(current_end, row_end)
        else:
            merged_intervals.append([current_chr, current_start, current_end])
            current_chr, current_start, current_end = row_chr, row_start, row_end

    merged_intervals.append([current_chr, current_start, current_end])
    return pd.DataFrame(merged_intervals, columns=['chr', 'start', 'end'])

def get_chip_contained_loci(chip_data, chr_, start, end):
    """ Checks if a ChIP-seq peak is fully contained within a loop locus and returns True/False. """
    if chip_data is None or chip_data.empty:
        return False  # Return False if there is no valid ChIP-seq data
    
    contained = chip_data[(chip_data['chr'] == chr_) & (chip_data['start'] >= start) & (chip_data['end'] <= end)]
    
    return not contained.empty  # Return True if at least one peak is found

# ==============================
# Process Each Model and Chromosome
# ==============================

summary_results = []
missing_files = 0
missing_chipseq_files = 0

for model in models:
    model_name = os.path.basename(model)  # Extract model name from path
    
    for chrom in chromosomes:
        print(f"Processing Model: {model_name}, Chromosome: {chrom}")
        # for enhanced one:
        loop_file = f"{model}/hiccups_results_{chrom}/merged_loops.bedpe"

        # # for LR one:
        # loop_file = f"{model}/hiccups_results_lr_{chrom}/merged_loops.bedpe"


        if not os.path.exists(loop_file):
            print(f"⚠ Warning: Skipping {loop_file}, file not found.")
            missing_files += 1
            continue  # Skip if loop file does not exist

        # Load loop callset
        try:
            HR_loop_callset = pd.read_csv(loop_file, sep='\t', skiprows=[1], usecols=[0, 1, 2, 3, 4, 5])
            HR_loop_callset.columns = ['#chr1', 'x1', 'x2', 'chr2', 'y1', 'y2']
        except Exception as e:
            print(f"❌ Error reading {loop_file}: {e}")
            continue

        print(f"✅ Loaded {loop_file}: {HR_loop_callset.shape[0]} loops")

        # Process Loop Loci (Expand by ±5kb)
        loop_callset_upstream = HR_loop_callset[['#chr1', 'x1', 'x2']].copy()
        loop_callset_downstream = HR_loop_callset[['chr2', 'y1', 'y2']].copy()
        loop_callset_downstream.columns = loop_callset_upstream.columns
        combined_loop_callset = pd.concat([loop_callset_upstream, loop_callset_downstream], ignore_index=True)
        combined_loop_callset = combined_loop_callset.drop_duplicates(keep='first')
        combined_loop_callset.columns = ['chr', 'start', 'end']
        combined_loop_callset = convert_chr_to_int(combined_loop_callset, loop_file)

        # Expand and merge intervals
        merged_loop_loci = expand_and_merge_intervals(combined_loop_callset)

        # Load ChIP-seq Data
        chipseq_dir = f"/GM12878_newsplit_chr18_22_test/TF_hg19_Cell/{cell_line}"
        chipseq_dfs = {}

        for factor in ["CTCF", "RAD21", "SMC3"]:
            chipseq_file = f"{chipseq_dir}/{factor}/merged_output.txt"
            if os.path.exists(chipseq_file) and os.path.getsize(chipseq_file) > 0:
                chipseq_dfs[factor] = pd.read_csv(chipseq_file, sep='\t', header=None, usecols=[0, 1, 2], names=['chr', 'start', 'end'])
                chipseq_dfs[factor] = convert_chr_to_int(chipseq_dfs[factor], chipseq_file)
                print(f"✅ Loaded {factor}: {chipseq_dfs[factor].shape[0]} peaks")
            else:
                print(f"⚠ Warning: {factor} data missing for {cell_line}")
                missing_chipseq_files += 1

        # Compute Overlap
        total_loci = len(merged_loop_loci)
        validated_loci = 0

        for _, row in merged_loop_loci.iterrows():
            chr_, start, end = row['chr'], row['start'], row['end']
            if pd.isna(chr_):
                continue  # Skip NaN chromosome rows

            available_factors = [factor for factor in chipseq_dfs if factor in chipseq_dfs]
            num_required_matches = len(available_factors)

            matched_count = sum(
                get_chip_contained_loci(chipseq_dfs[factor], chr_, start, end)
                for factor in available_factors
            )

            if matched_count == num_required_matches:
                validated_loci += 1

        validated_percentage = (validated_loci / total_loci) * 100 if total_loci > 0 else 0

        summary_results.append({
            'Cell_Line': cell_line,
            'Model': model_name,
            'Chromosome': chrom,
            'Total_Loci': total_loci,
            'Validated_Loci': validated_loci,
            'Validated_Percentage': round(validated_percentage, 2)
        })

# Save Final Summary
summary_df = pd.DataFrame(summary_results)
output_file = f"/Loop_{cell_line}_LR_validation_Feb11.csv"
summary_df.to_csv(output_file, index=False)

print(f"Final summary saved to: {output_file}")
