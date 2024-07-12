#This script reads in the table created for the phenotype we want to make PRSs for 
#It also loads the SNPs we have availble in FGFP so that we only make PRSs with these SNPs
#It outputs SNP files ready to be processed by PRSice

#Due to accessing the eur reference this file should be run on bp

#setwd to scripts/PhenoPRS
library(gwasrapidd)
library(dplyr)
library(ieugwasr)

#Data path for proj dir
data.path <- "../../../data/PhenoPRS"

#Also load in path to scratch in case need to access files here
source("../../parameters/scratch_path.txt")

#Lets load in the GWAS table we have generated
gwas_table <- read.csv(file = file.path(data.path, "PhenoGWASSum.csv"), header = T)

#Load in fgfp snps before loop for computational efficency
load(file.path(data.path, "all_snpswrsids.RData"))

#Start with the first GCST number to get the pipeline running
for (i in 1:nrow(gwas_table)){
  
  #Check if the code has already been run, otherwise can be skipped
  outputdir <- file.path(data.path, "prs_out", gwas_table$file_name[i])
  
  if(!dir.exists(outputdir)){
    
study <- get_studies(study_id = gwas_table$GCST[i])
snps <- get_associations(study_id = gwas_table$GCST[i])

#Make a SNPs df to be used for downstream analyses

snps_df <- data.frame(rsid = snps@risk_alleles[["variant_id"]],
                      eff.allele = snps@risk_alleles[["risk_allele"]],
                      EAF = snps@risk_alleles[["risk_frequency"]],
                      beta = snps@associations[["beta_number"]] * ifelse(snps@associations[["beta_direction"]] == "increase", 1, -1), #Multiply so directionality is correct
                      SE = snps@associations[["standard_error"]],
                      pval = snps@associations[["pvalue"]])

#Filter to just those snps we are interested in
merged <- filter(merged, rsid %in% snps_df$rs)

snps_df <- merge(snps_df, merged[c("V1", "rsid")], by = "rsid", all.y = T) %>% 
  rename(fgfp_linker = V1)

#Now we can clump for independent hits
clumped <- ld_clump(dplyr::tibble(rsid=snps_df$rsid,pval=snps_df$pval),
  plink_bin = genetics.binaRies::get_plink_binary(),
  bfile = file.path(scratch, "ldclump/EUR")
)

#Final dataframe with clumped SNPs
snps_df <- filter(snps_df, rsid %in% clumped$rsid)

#Make the final df and save output if the file does not already exist


if(!dir.exists(outputdir)){
  dir.create(outputdir)
}

outputfile <- file.path(outputdir, paste0(gwas_table$file_name[i], "_processedgwas.tsv"))

if(!file.exists(outputfile)){
  write.table(snps_df, file = outputfile, quote = F, row.names = F)
}


#Also write out a file with just the snps we are interested in for FGFP
outputfile_snps <- file.path(outputdir, paste0(gwas_table$file_name[i], "_snps.txt"))

if(!file.exists(outputfile_snps)){
  write.table(snps_df["fgfp_linker"], file = outputfile_snps, col.names=F, row.names=F, quote=F)
}

} else {cat(gwas_table$file_name[i], "has already been run \n")}

  }


