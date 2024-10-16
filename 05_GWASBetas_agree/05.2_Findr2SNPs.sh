#Bash script to find correlated variants with SNPs of interest
#This should be run on epifrank

# Save the current working directory
wd=$(pwd)

# Data directory and traits
data_dir=../../../data/PhenoPRS
traits=("triglycerides" "Serum_urea_levels" "Serum_uric_acid_levels" "yengo_bmi" "Gamma_glutamyltransferase_levels_UKB_data_field_30730_" "Alanine_aminotransferase_levels" "HDL_cholesterol" "Standing_height" "Hip_circumference_UKB_data_field_49_" "Red_blood_cell_count")

# Loop through each trait
for trait in "${traits[@]}"; do
    # Change to trait SNP directory
    cd "$data_dir/phenosumstats/$trait" || { echo "Directory not found for trait: $trait"; continue; }

    # Get rsids in new file
    tail -n +2 "${trait}_processedgwas.tsv" | awk '{print $1}' > rsid_list.txt

    # Remove duplicates
    sort -u rsid_list.txt -o rsid_list.txt

    # Now find correlated SNPs with our SNPs of interest

    # First make a new directory
    mkdir -p corsnps  # Use -p to avoid errors if the directory already exists

    plink --bfile ~/clump/EUR \
          --ld-snp-list rsid_list.txt \
          --ld-window-kb 500 \
          --ld-window 99999 \
          --ld-window-r2 0.01 \
          --r2 \
          --out "corsnps/corsnps_$trait"  # Include the trait name in the output file

    # Set back to original working directory
    cd "$wd" || exit  # Exit if returning to the original directory fails
done
