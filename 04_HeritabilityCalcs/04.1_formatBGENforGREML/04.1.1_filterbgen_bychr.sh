#!/bin/bash
#SBATCH --job-name="bgen_qc"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1 
#SBATCH --time=01:00:00
#SBATCH --mem=20G
#SBATCH --account=sscm009461
#SBATCH --array=1-22

################################################
#Script to calculate heritability with GREML
################################################

#Add in modules needed
module add lang/r/4.3.0-bioconductor-gcc
module add apps/qctool/2.0.7


#Extract slurmtaskid for chr number we are looking at
id=$SLURM_ARRAY_TASK_ID

# Add a leading zero for ids 1-9
if [ $id -lt 10 ]; then
  chr="0$id"
else
  chr="$id"
fi

# Echo to verify the chr value
echo "Editing chromosome: $chr"

#################################################
# Step 1: Filter SNPs to use for h2 calc
#################################################

#Extract our basedir for bp if needed
#source ../../../parameters/scratch_path.txt

#Assign data loc object
data_dir="../../../../data"

## First we identify SNPs with INFO >= 0.9 & MAF > 0.05 and HWE > 1e-5
## From the previous PRS generation project we have the list of SNPs that are in the imputed bgen file, this will be used for SNP identification, filtering on the parameters above

#SNP_stat files are located here:
snpstat_dir="${data_dir}/FGFP_data/genetic_data/prsfgfp/bgen_w_rsids/snp_stat"
#This is where the data will be go
output_dir="${data_dir}/PhenoPRS/greml/bgen"


#Make directory for output of Rscript filter
mkdir $output_dir/snps2use

#Make input and output for Rscript
input="$snpstat_dir/snp_summary${chr}.txt"
output_snps="$output_dir/snps2use/QC_snps${chr}.txt"
#Call Rscript to filter our specified values and save 
Rscript 04.1.2_filtersnpstat.R $input $output_snps

#####################################################
# Step 2: Use SNPs file generated to filter bgen file
#####################################################
input_bgen="${data_dir}/FGFP_data/genetic_data/prsfgfp/bgen_w_rsids/bgen_w_rsids/data_chr${chr}.bgen"
output_bgen="${data_dir}/PhenoPRS/greml/bgen/qc_chr${chr}.bgen"


qctool -g "$input_bgen" -og "$output_bgen" -incl-rsids "$output_snps"


################################################################
# Step 3: Add Chr number to this file as it is currently missing
#################################################################
output_bgen_wchr="${data_dir}/PhenoPRS/greml/bgen/qc_chr${chr}_wchr.bgen"

qctool -g "$output_bgen" -og "$output_bgen_wchr" -assume-chromosome ${chr}

#Now we can remove the previous bgen
rm $output_bgen


