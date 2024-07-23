#### Now we can run bivariate GREML with our .phen file and GRM

#Set working directory to where the script is located
#Set a base_dir
geno_prefix="../../../../data/PhenoPRS/greml"

#Add greml
module add apps/gcta/1.93.2-beta

#First extract the pheno file to extract phenotypes and number of loops need to run
pheno_file="${geno_prefix}/phenos4greml/bivariateGREML/phen4cor.phen"

# Get the headers of the columns in the phenotypes file
pheno_names=${geno_prefix}/phenos4greml/bivariateGREML/phennames.txt

# Determine the number of columns
num_cols=$(awk -F'\t' '{print NF}' $pheno_file | head -n 1)


#Lets try to compute h2 with our rm_025 and grm file
# Loop through values 1 to 7 as this is the number of phenotypes in the file
for num in $(seq 3 2 $num_cols); do
  echo "Running GCTA for mpheno = $num"

  index=$((num - 2))
  index2=$((index + 1))
  pheno1=$(sed -n "${index}p" "$pheno_names")
  pheno2=$(sed -n "${index2}p" "$pheno_names")

  pheno1__pheno2=${pheno1}__${pheno2}
  
  gcta64 --grm ${geno_prefix}/greml_out/grm_rm025/grm_rm025 \
         --pheno $pheno_file \
         --reml-bivar $index $index2 \
         --out ${geno_prefix}/greml_out/h2out_cor/h2out_$pheno1__pheno2 \
         --thread-num 10
done

