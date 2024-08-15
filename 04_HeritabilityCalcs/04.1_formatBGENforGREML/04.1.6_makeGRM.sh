### This script calculates the grm from our genetic data


#Set working directory to where the script is located

#Set a base_dir
geno_prefix="../../../../data/PhenoPRS/greml"

#Add greml
module add apps/gcta/1.93.2-beta
module add apps/plink/2.0.0


#Make the GRM with no removal of people
gcta64 --bfile ${geno_prefix}/plink/pruned/pruned_data --make-grm --out ${geno_prefix}/greml_out/grm/grm --thread-num 10

#This code removes cryptic relatedness if needs:
gcta64 --grm ${geno_prefix}/greml_out/grm/grm --grm-cutoff 0.025 --make-grm --out ${geno_prefix}/greml_out/grm_rm025/grm_rm025


#### Now we can make a grm with daves sample size
# First use plink to remove the ids we are not interested in

plink2 --bfile ${geno_prefix}/plink/pruned/pruned_data --remove ${geno_prefix}/plink/pruned/dave_sample/relatedness_ids2remove.txt --make-bed --out ${geno_prefix}/plink/pruned/dave_sample/dave_sample
gcta64 --bfile ${geno_prefix}/plink/pruned/dave_sample/dave_sample --make-grm --out ${geno_prefix}/greml_out/grm/grm_davesample --thread-num 10

