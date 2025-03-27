#!/bin/bash

# List of chromosome numbers
CHR_NUMS=(16 17 18 19)

# List of folders to process (update these paths as needed)
FOLDERS=(
    "/your/path/to/results"
    # Add more folders here as needed
)

# Loop through each folder
for folder in "${FOLDERS[@]}"; do
    echo "Processing folder: $folder"

    # Change to the directory of the current folder
    cd "$folder" || { echo "Cannot access folder $folder"; continue; }

    # Loop through each chromosome number
    for CHR_NUM in "${CHR_NUMS[@]}"; do
        echo "Processing chromosome ${CHR_NUM} in folder ${folder}..."

        # Run the HiCCUPs command for LR data
        java -jar juicer_tools_1.22.01.jar hiccups -r 10000 -k KR \
        total_merged_10kb.hic \
        -c chr${CHR_NUM} \
        hiccups_results_ori_KR_chr${CHR_NUM}
    done

    echo "Completed processing for folder: $folder"
done

echo "All folders and chromosomes processed."