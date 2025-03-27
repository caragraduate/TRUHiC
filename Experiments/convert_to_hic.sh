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

    # Loop through each chromosome number
    for CHR_NUM in "${CHR_NUMS[@]}"; do
        echo "Processing chromosome ${CHR_NUM} in folder ${folder}..."

        # Modify the input and output filenames for each chromosome of human
        java -Xmx1024g -jar juicer_tools_1.22.01.jar pre -r 10000 \
        preds_lr_test_chr${CHR_NUM}_ratio${ratio}_convert.txt \
        preds_lr_test_chr${CHR_NUM}_ratio${ratio}_convert.hic hg19

        # # Modify the input and output filenames for each chromosome of mouse
        # java -Xmx1024g -jar juicer_tools_1.22.01.jar pre -r 10000 \
        # preds_lr_test_chr${CHR_NUM}_ratio16_convert.txt \
        # preds_lr_test_chr${CHR_NUM}_ratio16_convert.hic mm9

        echo "Finished .hic conversion for chromosome ${CHR_NUM} in folder ${folder}."

        # Run the arrowhead command
        java -Xmx1024g -jar juicer_tools_1.22.01.jar arrowhead -r 10000 -k KR \
        preds_lr_test_chr${CHR_NUM}_ratio${ratio}_convert.hic \
        preds_lr_test_chr${CHR_NUM}_ratio${ratio}_convert_10kb --threads 16

        echo "Finished arrowhead processing for chromosome ${CHR_NUM} in folder ${folder}."
    done

    echo "Completed processing for folder: $folder"
done

echo "All folders and chromosomes processed."