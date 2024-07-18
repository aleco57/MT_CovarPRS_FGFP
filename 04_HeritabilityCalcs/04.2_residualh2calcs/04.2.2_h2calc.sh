
cd /user/work/am17168/greml_calc

module add apps/gcta/1.93.2-beta

base_dir="/user/work/am17168/greml_calc"

#Make the GRM with no removal of people
gcta64 --bfile ${base_dir}/data/plink/pruned/pruned_data --make-grm --out ${base_dir}/data/greml_out/grm/grm --thread-num 10

#This code removes cryptic relatedness if needs:
gcta64 --grm ${base_dir}/data/greml_out/grm/grm --grm-cutoff 0.025 --make-grm --out ${base_dir}/data/greml_out/grm_rm025/grm_rm025

#Lets try to compute h2 with our rm_025 and grm file
# Loop through values 1 to 7 as this is the number of phenotypes in the file
for num in {1..7}; do
  echo "Running GCTA for mpheno = $num"
  
  gcta64 --grm ${base_dir}/data/greml_out/grm_rm025/grm_rm025 \
         --pheno ${base_dir}/data/mapping_files/phenos4greml/sigRNT_plinkids.phen \
         --mpheno $num \
         --reml \
         --out ${base_dir}/data/greml_out/h2out/h2out_pheno$num \
         --thread-num 10
done


#Now we can do the same for our PRS adjusted phenotypes
# Loop through values 1 to 13
for num in {1..13}; do
  echo "Running GCTA for mpheno = $num"
  
  gcta64 --grm ${base_dir}/data/greml_out/grm_rm025/grm_rm025 \
         --pheno ${base_dir}/data/mapping_files/phenos4greml/sigRNT_PRSRes_plinkids.phen \
         --mpheno $num \
         --reml \
         --out ${base_dir}/data/greml_out/h2out/prsadj/h2outadj_pheno$num \
         --thread-num 10
done


#Now copy these heritability files over to proj 007 for processing there


