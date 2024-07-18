#!/bin/bash
#SBATCH --job-name="plink_convert"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1 
#SBATCH --time=03:00:00
#SBATCH --mem=40G
#SBATCH --account=sscm009461

#Add in modules needed
#module add lang/r/4.3.0-bioconductor-gcc
module add apps/qctool/2.0.7
module add apps/plink/2.0.0

#Set a base_dir
geno_prefix="../../../../data/PhenoPRS/greml"


##################################################################################
# Convert to plink for hard-thresh probability coding
##################################################################################

#Make directory for plink output
mkdir ${geno_prefix}/plink

#Convert to plink file format
qctool -g ${geno_prefix}/bgen/qcwchr_allchr.bgen -ofiletype binary_ped -og ${geno_prefix}/plink/data_allchr

##################################################################################
# LD prune using plink
##################################################################################
plink2 --bfile ${geno_prefix}/plink/data_allchr --indep-pairwise 50 5 0.45 --out ${geno_prefix}/plink/pruned/pruned_data

# Step 2: Create a new dataset with pruned SNPs
plink2 --bfile ${geno_prefix}/plink/data_allchr --extract ${geno_prefix}/plink/pruned/pruned_data.prune.in --make-bed --out ${geno_prefix}/plink/pruned/pruned_data
