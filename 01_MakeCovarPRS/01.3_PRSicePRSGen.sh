#!/bin/bash
#Note this script must be run in bp
#Set wd to /scripts/MT_CovarPRS_FGFP/01_MakeCovarPRS
wd=$(pwd)

#Extract our basedir for bp
source ../../parameters/scratch_path.txt

# Use find to get all directories in the target directory and loop over them
# Set the path to your text file
file_with_directories="../../../data/PhenoPRS/phenosumstats/phenosumstats_dirnames.txt"

# Loop through each line (directory name) in the file
while IFS= read -r dir_name; do
        # Process each directory name
        echo "Processing directory: $dir_name"
        
       cd ../../../data/PhenoPRS/phenosumstats/${dir_name}

    #Make output dir
    mkdir prsice_out

    #Make fgfp_basedir
    fgfpdir="../../../FGFP_data"

#Note a random outcome is chosen just for instrument to be developed
    Rscript ~/PRSice/PRSice.R \
        --prsice ~/PRSice/./bin/PRSice \
        --base ${dir_name}_processedgwas.tsv \
        --snp fgfp_linker \
        --A1 eff.allele  \
        --pvalue pval \
        --stat beta \
        --beta \
        --target ${fgfpdir}/genetic_data/prsfgfp/qc_bgen/withchr/qc_bgen# \
        --type bgen \
        --pheno-col class.Actinobacteria.id.419 \
        --pheno ${fgfpdir}/FGFP_MT/mt_covadj_allids.txt \
        --thread 10 \
        --binary-target F \
        --ignore-fid \
        --all-score \
        --quantile 10 \
        --out prsice_out/${dir_name} \
        --allow-inter \
        --bar-levels 1 \
        --no-clump \
        --extract ${dir_name}_snps.txt

    cd prsice_out/
    rm *.png

#Set back to wd so can loop back around
cd $wd

done < "$file_with_directories"
    
 


   