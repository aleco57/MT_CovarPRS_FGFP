#### Now we can estimate h2 with our .phen file and GRM

#Set working directory to where the script is located
#Set a base_dir
geno_prefix="../../../../data/PhenoPRS/greml"

#Add greml
module add apps/gcta/1.93.2-beta

#First extract the pheno file to extract phenotypes and number of loops need to run
pheno_file="${geno_prefix}/phenos4greml/h2mts.phen"

# Get the headers of the columns in the phenotypes file
pheno_names=${geno_prefix}/phenos4greml/phennames.txt

# Determine the number of columns
num_cols=$(awk -F'\t' '{print NF}' $pheno_file | head -n 1)


#Lets try to compute h2 with our rm_025 and grm file
# Loop through values 1 to 7 as this is the number of phenotypes in the file
for num in $(seq 3 $num_cols); do
  echo "Running GCTA for mpheno = $num"

  index=$((num - 2))
  pheno=$(sed -n "${index}p" "$pheno_names")
  
  
  gcta64 --grm ${geno_prefix}/greml_out/grm_rm025/grm_rm025 \
         --pheno ${geno_prefix}/phenos4greml/h2mts.phen \
         --mpheno $index \
         --reml \
         --out ${geno_prefix}/greml_out/h2out/h2out_$pheno \
         --thread-num 10
done

