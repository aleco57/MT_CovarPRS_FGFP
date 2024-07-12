#!/bin/bash
#Note this script must be run in bp
#Set wd to /scripts/MT_CovarPRS_FGFP/01_MakeCovarPRS
wd=$(pwd)

#Extract our basedir for bc
source ../../parameters/scratch_path.txt

#Extract the phenos for us to loop round
#Need to update this to a loop !
num_rows=$(awk 'END{print NR}' ../../../data/PhenoPRS/PhenoGWASSum.csv)

# Loop from NR==2 to the number of rows in the CSV file
for ((i=2; i<=num_rows; i++)); do
    trait=$(awk -F ',' -v row=$i 'NR==row {print $3}' ../../../data/PhenoPRS/PhenoGWASSum.csv)
    # Use 'trait' variable here for further processing
    echo "Trait for row $i: $trait"

trait_dir="../../data/PhenoPRS/${trait}"

if [ ! -d "$trait_dir" ]; then
    echo "Please run script 01 for trait: ${trait}"
elif [ -d "$trait_dir/prsice_out" ]; then
    echo "PRSice has already been run for trait: ${trait}"
else
    cd ../../../data/PhenoPRS/prs_out/${trait}

    #Make output dir
    mkdir prsice_out

    #Make fgfp_basedir
    fgfpdir="../../../data/FGFP_data"

#Note a random outcome is chosen just for instrument to be developed
    Rscript ~/PRSice/PRSice.R \
        --prsice ~/PRSice/./bin/PRSice \
        --base ${trait}_processedgwas.tsv \
        --snp fgfp_linker \
        --A1 eff.allele  \
        --pvalue pval \
        --stat beta \
        --beta \
        --target ${fgfpdir}/genetic_data/prsfgfp/qc_bgen/withchr/qc_bgen# \
        --pheno-col class.Actinobacteria.id.419 \
        --type bgen \
        --pheno ${fgfpdir}/FGFP_MT/mt_covadj_allids.txt \
        --thread 10 \
        --binary-target F \
        --ignore-fid \
        --all-score \
        --quantile 10 \
        --out prsice_out/${trait} \
        --allow-inter \
        --bar-levels 1 \
        --no-clump \
        --extract ${trait}_snps.txt

    cd prsice_out/
    rm *.png

    echo "PRSice run for trait: ${trait}"

fi

cd $wd

done
