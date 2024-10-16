#### This script will look at our candidate relationships and see if these covariate SNPs are showing in the sumstats for those specific traits in mibiogen

#Set up
#Library
library(tidyverse)
library(data.table)
library(qqman)
library(ggrepel)
library(LDlinkR)
library(broom)


#Data path for proj dir
data.path <- "../../../data/PhenoPRS"

#Also load in path to scratch, proj_dir and fgfp data
source("../../parameters/base_dir.R")
source("../../parameters/scratch_path.txt")

#Load in the bus we are interested in 
load(file.path(data.path, "data_out/data_ConsistentAssoc03.2.RData"))

#Also load in the SNP data
load(file.path(data.path, "phenosumstats/filtered_sumstats_list.RData"))

load(file.path(data.path, "data_out/filtered_bug_sumstats_05.1.RData"))


#Lets look at FGFP sumstats
############################

#Get our candidiates we are intersted in
#
lm_out <- c()
for(covar in names(merged_gwasbetas)){
  for(bug in names(merged_gwasbetas[[covar]])){
    out <- lm(frequentist_add_beta_1 ~ beta, data = merged_gwasbetas[[covar]][[bug]]) %>% summary() %>% tidy() %>% filter(term == "beta")
    
    out <- cbind(out, covar, bug)
    lm_out <- rbind(lm_out, out)
  }}

#Lets focus on those bugs where we saw the replication in beta coefficients but could look at all our candidiate relationships
candidates <- filter(lm_out, p.value < 0.05/10)

#Lets also just focus on taxa traits, we can look at alpha div with MiBioGen
candidates <- candidates[!grepl("Div_", candidates$bug),]

#Rename the names of snps_df for consistent naming
#And specify the function for generating the file names
#Specify a function for cleaning directory names
clean_directory_name <- function(x) {
  x %>%
    stringr::str_replace_all("[^[:alnum:]]+", "_") %>% 
    stringr::str_replace_all("_+", "_")   
}

names(snps_df) <- clean_directory_name(names(snps_df))



manhattan_results <- list()

for(i in 1:nrow(lm_out)){
  print(lm_out[i,])
  
  covar <- lm_out$covar[i]
  bug <- lm_out$bug[i]
  
  snpsOfInterest <- snps_df[[covar]]$rsid
  
  #Read in the bug sumstat
  bug <- gsub("Res", "", bug)
  bug_sumstats <- data.table::fread(file = file.path(data.path, "davesRNTsumstats", paste0(bug, "_allchr.txt.gz")))
  
  
  #Filter to only those SNPs which are present in the sumstat
  snpsOfInterest <- snpsOfInterest[snpsOfInterest %in% bug_sumstats$rsid]
  
  #Filter out bad SNPs
  bug_sumstats_clean <- bug_sumstats[!is.na(bug_sumstats$frequentist_add_pvalue), ]
  
  #Make output file
  output_file <- file.path(data.path, "data_out/manhattan/manhattan_fgfp", paste0(covar, "_", bug, "_manhattan_plot1.png"))
  
  # Open a PNG device
  png(filename = output_file, width = 1200, height = 800, res = 150)
  
  #Now we can plot our manhattan plot
  manhattan(bug_sumstats_clean, 
            chr="chromosome", 
            bp="position", 
            snp="rsid", 
            p="frequentist_add_pvalue", 
            highlight = snpsOfInterest)
  
  
  # Close the PNG device
  dev.off()
  
  #Now we can plot also those SNPs which are in LD with our SNPs of interest
  #First read in the file
  ldsnps_file <- file.path(data.path, "phenosumstats", covar, "corsnps", paste0("corsnps_", covar, ".ld"))
  
  ldsnps <- read.table(file= ldsnps_file, stringsAsFactors = FALSE, header = TRUE)
  
  snpsOfInterest2 <- ldsnps$SNP_B %>% unique()
  snpsOfInterest2 <- snpsOfInterest2[snpsOfInterest2 %in% bug_sumstats$rsid]
  
  #Plot
  
  #Make output file
  output_file2 <- file.path(data.path, "data_out/manhattan/manhattan_fgfp", paste0(covar, "_", bug, "_manhattan_plot2.png"))
  
  # Open a PNG device
  png(filename = output_file2, width = 1200, height = 800, res = 150)
  
  
  manhattan(bug_sumstats_clean, 
            chr="chromosome", 
            bp="position", 
            snp="rsid", 
            p="frequentist_add_pvalue", 
            highlight = snpsOfInterest2)
  
  # Close the PNG device
  dev.off()
  
  #Lets also extract those SNPs which reach the relaxed threshold and see if the betas are correlated in the direction we expect
  manhattan_results[[paste0(covar, "_", bug)]][["sigsnps"]] <- filter(bug_sumstats, rsid %in% snpsOfInterest2, frequentist_add_pvalue < 1e-5)
}


#Save the results so can visualise in Rstudio
save(manhattan_results, file = file.path(data.path, "data_out/manhattan/manhattan_fgfp/sig_snps_05.3.RData"))



####### Lets also look at the BMI snps that were genome-wide sig for the G__unclassified_P_Firmicutes_RNTRes bug
manhattan_results[["yengo_bmi_G_unclassified_P_Firmicutes_RNT"]]$sigsnps
covar <- "yengo_bmi"
ldsnps_file <- file.path(data.path, "phenosumstats", covar, "corsnps", paste0("corsnps_", covar, ".ld"))

ldsnps <- read.table(file= ldsnps_file, stringsAsFactors = FALSE, header = TRUE)

filter(ldsnps, SNP_B %in% manhattan_results[["yengo_bmi_G_unclassified_P_Firmicutes_RNT"]]$sigsnps$rsid)[["R2"]]

highldsnps <- filter(ldsnps, SNP_B %in% manhattan_results[["yengo_bmi_G_unclassified_P_Firmicutes_RNT"]]$sigsnps$rsid & R2 > 0.1)[["SNP_B"]]