import os
import sys
import math
import numpy as np
import pandas as pd

# List of subfolder paths to process
data_folders = [
    "/your/path/to/results",
    # Add more folders here as needed
]

# List of chromosomes to process
chromosomes = ['chr18', 'chr19', 'chr20', 'chr21', 'chr22']
ratio = 16  # Ratio to use (update as needed)
res = 10000  # Resolution
size = 40  # Size of the crop for Hi-C matrix

# Print all data folders to be processed
print("List of data folders to be processed:")
for folder in data_folders:
    print(f"  - {folder}")

def crop_hic_matrix_by_chrom(chrom, hr_contacts_dict, lr_contacts_dict, size=40, thred=200):
    distance = []
    crop_mats_hr = []
    crop_mats_lr = []
    indices = []
    row, col = hr_contacts_dict[chrom].shape

    if row <= thred or col <= thred:
        print('HiC matrix size wrong!')
        sys.exit()

    def quality_control(mat, thred=0.05):
        return len(mat.nonzero()[0]) >= thred * mat.shape[0] * mat.shape[1]

    for idx1 in range(0, row - size, size):
        for idx2 in range(0, col - size, size):
            if abs(idx1 - idx2) < thred:
                if quality_control(lr_contacts_dict[chrom][idx1:idx1+size, idx2:idx2+size]):
                    distance.append([idx1 - idx2, chrom])

                    lr_contact = lr_contacts_dict[chrom][idx1:idx1+size, idx2:idx2+size]
                    hr_contact = hr_contacts_dict[chrom][idx1:idx1+size, idx2:idx2+size]

                    crop_mats_lr.append(lr_contact)
                    crop_mats_hr.append(hr_contact)
                    indices.append((idx1, idx2))

    crop_mats_hr = np.concatenate([item[np.newaxis, :] for item in crop_mats_hr], axis=0)
    crop_mats_lr = np.concatenate([item[np.newaxis, :] for item in crop_mats_lr], axis=0)

    return crop_mats_hr, crop_mats_lr, distance, indices

def write_matrices_to_txt(hr_mats, indices, filename, res=10000):
    entries = {}  # Dictionary to store entries as (pos1, pos2): [score_sum, count]

    for i in range(len(hr_mats)):
        idx1, idx2 = indices[i]
        for row_idx in range(hr_mats.shape[1]):  # Iterate over rows
            for col_idx in range(hr_mats.shape[2]):  # Iterate over columns
                score = hr_mats[i, row_idx, col_idx, 0]  # Get the score

                if score != 0:  # Exclude lines with a score of zero
                    pos1 = (idx1 + row_idx) * res
                    pos2 = (idx2 + col_idx) * res

                    # Use sorted order to ensure (pos1, pos2) and (pos2, pos1) are the same
                    key = tuple(sorted((pos1, pos2)))

                    # Accumulate scores and counts for averaging
                    if key in entries:
                        entries[key][0] += score  # Sum the scores
                        entries[key][1] += 1  # Count the occurrences
                    else:
                        entries[key] = [score, 1]  # Initialize with the first score and count of 1

    # Compute the average score for each pair
    averaged_entries = [(key[0], key[1], score_sum / count) for key, (score_sum, count) in entries.items()]

    # Sort entries based on pos2 and pos1
    sorted_entries = sorted(averaged_entries, key=lambda x: (x[1], x[0]))

    # Write sorted entries to file
    with open(filename, 'w') as f:
        for entry in sorted_entries:
            f.write(f"{entry[0]}\t{entry[1]}\t{entry[2]}\n")

# Loop through all subfolders and chromosomes
for data_file in data_folders:
    print(f"\nProcessing folder: {data_file}")
    
    for chrom in chromosomes:
        print(f"\nProcessing chromosome: {chrom}")
        enhanced_hic = os.path.join(data_file, f'preds_lr_test_{chrom}_ratio{ratio}.npy')

        # Print which file is currently being processed
        print(f"Loading enhanced Hi-C data from: {enhanced_hic}")
        our_test = np.load(enhanced_hic)

        # Chromosome length details
        chrom_len_file = 'chromosome.txt'
        # chrom_len_file = '/mm9.chrom.sizes.txt'

        chrom_len = {item.split()[0]: int(item.strip().split()[1]) for item in open(chrom_len_file).readlines()}
        mat_dim = int(math.ceil(chrom_len[chrom] * 1.0 / res))

        # Process hr_hic and lr_hic matrices
        GM12878_file= '/GM12878/intra_NONE'
        GM12878_rep_file= '/GM12878_4/intra_NONE'

        hr_hic_file = os.path.join(GM12878_file, f'{chrom}_10k_intra_NONE.txt')
        lr_hic_file = os.path.join(GM12878_rep_file, f'{chrom}_10k_intra_NONE_downsample_ratio{ratio}.txt')

        print(f"Processing high-resolution Hi-C file: {hr_hic_file}")
        hr_contacts_dict = {}
        hr_contact_matrix = np.zeros((mat_dim, mat_dim))
        for line in open(hr_hic_file).readlines():
            idx1, idx2, value = map(float, line.strip().split('\t')[:3])
            if idx2 / res >= mat_dim or idx1 / res >= mat_dim:
                continue
            hr_contact_matrix[int(idx1 / res)][int(idx2 / res)] = value
        hr_contact_matrix += hr_contact_matrix.T - np.diag(hr_contact_matrix.diagonal())
        hr_contacts_dict[chrom] = hr_contact_matrix

        print(f"Processing low-resolution Hi-C file: {lr_hic_file}")
        lr_contacts_dict = {}
        lr_contact_matrix = np.zeros((mat_dim, mat_dim))
        for line in open(lr_hic_file).readlines():
            idx1, idx2, value = map(float, line.strip().split('\t')[:3])
            if idx2 / res >= mat_dim or idx1 / res >= mat_dim:
                continue
            lr_contact_matrix[int(idx1 / res)][int(idx2 / res)] = value
        lr_contact_matrix += lr_contact_matrix.T - np.diag(lr_contact_matrix.diagonal())
        lr_contacts_dict[chrom] = lr_contact_matrix

        # Crop matrices
        crop_mats_hr, crop_mats_lr, distance, indices = crop_hic_matrix_by_chrom(chrom, hr_contacts_dict, lr_contacts_dict)

        # Write matrices to text file
        enhanced_hic_txt = os.path.join(data_file, f'preds_lr_test_{chrom}_ratio{ratio}.txt')
        print(f"Writing matrices to text file: {enhanced_hic_txt}")
        write_matrices_to_txt(our_test, indices, enhanced_hic_txt)

        # Read the existing DataFrame
        output_matrices_chr = pd.read_csv(enhanced_hic_txt, names=['pos1', 'pos2', 'score'], sep='\t')
        output_matrices_chr.loc[11]

        chromosome = enhanced_hic_txt.split('/')[-1].split('_')[-2].split('.')[0]
        # for hg19, we don't need the "chr" prefix so we need the code in the next line
        # chromosome = chromosome[3:]
        ## for mouse, we DO need the "chr" prefix, so we don't need the command above, we need to comment # the line above
        print(chromosome)
        
        # Format the data
        generated_data = [f'0 {chromosome} {pos1} 0 0 {chromosome} {pos2} 1 {score}' 
                        for pos1, pos2, score in zip(output_matrices_chr['pos1'], output_matrices_chr['pos2'], output_matrices_chr['score'])]

        print(len(generated_data))

        enhanced_hic_txt_convert = data_file+ '/preds_lr_test_%s_ratio%d_convert.txt'%(chrom,ratio)

        with open(enhanced_hic_txt_convert, 'w') as f:
            for line in generated_data:
                f.write(line + '\n')

    print(f"Completed processing for folder: {data_file}\n")

print("All folders and chromosomes processed.")