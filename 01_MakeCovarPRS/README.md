## 01_MakeCovarPRS Contents

This folder contains the first step in the pipeline where microbiome PRSs are developed
In this folder there are four scripts and one subdirectory, there are:

(01) This is the .sh file to submit the 01.1 and 01.2 scripts to SLURM
(01.1) This is an R script which extracts data from a csv file (manually created from Fallony study) to extract GWAS data from the GWAS Catalog and filter to those SNPs we have availble in FGFP
(01.2) This is a shell script which uses FGFP genetic data in bgen file format and GWAS meta-data to generate the PRSs for the covariables in FGFP using the PRSice-2 wrapper
(01.3) Some issues when generating PRSs using the above scripts and GWAS Catalog, therefore subsequent PRSs were made using OpenGWAS as were set up better for PRS gen. Also added a function in this script to calculate the number of SNPs associated to our covariates before and after clumping

(ManualPRSGen) This is a subdirectory for making one of the PRSs manually by converting to gen file and generating a bespoke R script. It was used as sanity check to make sure PRSs were being generated correctly with the PRSice algorithm


### Known issues / Check
1. Slight differences in PRS when manually generated. There are almost identical however issues could be due to PRSice not converting to dosages properly of missingness accounted for differently when calculating manually (MAF is imputed for missing data in PRSice). Would be good to check this when get a chance 
