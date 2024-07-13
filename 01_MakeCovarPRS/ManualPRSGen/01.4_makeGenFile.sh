### For our creatinine variable we are seeing an inverse association with our trait of interest, these will be manually generated to see what is happening
### It will also confirm if we get compariable results with PRSice-2 algorithm 

#####################################################################
# Use QCtools to make gen files for our traits of interest
#####################################################################
#locate to script wd

module add apps/qctool/2.0.7

### Extract the snps of interest and generate gen files using QC tools
trait="creatinine3"
data_dir="../../../../data/PhenoPRS"
genetic_dir="../../../../data/FGFP_data/genetic_data/prsfgfp"

bgen_file="${genetic_dir}/qc_bgen/genotype_file.bgen"
sample_file="${genetic_dir}/pheno/sampleids.sample" 
snplist_file="${data_dir}/prs_out/${trait}/${trait}_snps.txt"

outfile="${data_dir}/prs_out/${trait}/${trait}"

qctool -g ${bgen_file} -s ${sample_file} -incl-rsids "${snplist_file}" -og ${outfile}.gen -os ${outfile}.sample 

#Now we have the .gen file, we can compute the PRS ourselves in R
#Note: the .gen file is structured with 3 numbers for each individual, representing the probability of AA,AB and BB alleles
#Which allele is A and which is B is determined by the ordering of the alleles on the two columns on .gen file
#We will calculate the PRS using the same methodology as PRSice2, which uses: beta*dosage/nAlleles
#Dosage(B) = 2*P(BB) + 1*P(AB)

#Rscript is located as file script 01.5

