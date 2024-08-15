#### Now we can estimate h2 with our .phen file and GRM

#Set working directory to where the script is located
#Set a base_dir
geno_prefix="../../../../data/PhenoPRS/greml"

#Add greml
module add apps/gcta/1.93.2-beta


###############################################
# Dircons phenotypes
###############################################

#First extract the pheno file to extract phenotypes and number of loops need to run
pheno_file="${geno_prefix}/phenos4greml/dircons/dircons.phen"

# Get the headers of the columns in the phenotypes file
pheno_names=${geno_prefix}/phenos4greml/dircons/dircons_phenos.txt

# Determine the number of columns
num_cols=$(awk -F'\t' '{print NF}' $pheno_file | head -n 1)


#Lets try to compute h2 with our rm_025 and grm file
# Loop through values 1 to 7 as this is the number of phenotypes in the file
for num in $(seq 3 $num_cols); do
  echo "Running GCTA for mpheno = $num"

  index=$((num - 2))
  pheno=$(sed -n "${index}p" "$pheno_names")
  
  
  gcta64 --grm ${geno_prefix}/greml_out/grm/grm_davesample \
         --pheno ${pheno_file} \
         --mpheno $index \
         --reml \
         --out ${geno_prefix}/greml_out/h2_out/dircons/h2out_$pheno \
         --thread-num 10
done


###############################################
# h2 bugs phenotypes
###############################################

#First extract the pheno file to extract phenotypes and number of loops need to run
pheno_file="${geno_prefix}/phenos4greml/h2bugs/h2bugs.phen"

# Get the headers of the columns in the phenotypes file
pheno_names=${geno_prefix}/phenos4greml/h2bugs/h2bugs.txt

# Determine the number of columns
num_cols=$(awk -F'\t' '{print NF}' $pheno_file | head -n 1)


#Lets try to compute h2 with our rm_025 and grm file
# Loop through values 1 to 7 as this is the number of phenotypes in the file
for num in $(seq 3 $num_cols); do
  echo "Running GCTA for mpheno = $num"

  index=$((num - 2))
  pheno=$(sed -n "${index}p" "$pheno_names")
  
  
  gcta64 --grm ${geno_prefix}/greml_out/grm/grm_davesample \
         --pheno ${pheno_file} \
         --mpheno $index \
         --reml \
         --out ${geno_prefix}/greml_out/h2_out/h2bugs/h2out_$pheno \
         --thread-num 10
done


