#!/bin/bash

# List of chromosome numbers
CHR_NUMS=(18 19 20 21 22)
ratio=16  # Ratio to use (update as needed)

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

    # Create the "fanc_insulation_KR" directory if it doesn't exist
    mkdir -p fanc_insulation_KR

    # Loop through each chromosome number
    for CHR_NUM in "${CHR_NUMS[@]}"; do
        echo "Processing chromosome ${CHR_NUM} in folder ${folder}..."

        # Check if the insulation file exists
        insulation_file="fanc_insulation_KR/inter30_10kb_insulation_chr${CHR_NUM}"

        if [ -f "$insulation_file" ]; then
            echo "Insulation file for chromosome ${CHR_NUM} already exists, skipping insulation calculation."
        else

            # Run the fanc insulation command
            fanc insulation preds_lr_test_chr${CHR_NUM}_ratio${ratio}_convert.hic@10kb@KR \
            fanc_insulation_KR/inter30_10kb_insulation_chr${CHR_NUM} \
            -w 100000 \
            -o bed

            echo "Finished insulation calculation for chromosome ${CHR_NUM} in folder ${folder}."
        
        fi
        
        # Run the fanc boundaries command
        fanc boundaries fanc_insulation_KR/inter30_10kb_insulation_chr${CHR_NUM} \
        fanc_insulation_KR/inter30_10kb_boundaries_100kb_chr${CHR_NUM} \
        -w 100kb -s 0.2

        echo "Finished TAD boundary processing for chromosome ${CHR_NUM} in folder ${folder}."
    done

    echo "Completed processing for folder: $folder"
done

echo "All folders and chromosomes processed."
