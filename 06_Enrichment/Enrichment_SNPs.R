### Clean script for finding covariate SNPs in microbiome GWAS results


#Lib
library("dplyr")
library("data.table")
library("RColorBrewer")
library("purrr")
library("TwoSampleMR")
library("ggplot2")
library(stringr)
library(purrr)
library(TwoSampleMR)

#Data path for proj dir
data.path <- "../../../data/PhenoPRS"

#Also load in path to scratch, proj_dir and fgfp data
source("../../parameters/base_dir.R")
source("../../parameters/scratch_path.txt")

#Also load in the SNP data
load(file.path(data.path, "phenosumstats/filtered_sumstats_list.RData"))

#Load in the candidate bugs
load(file.path(data.path, "data_out/data_ConsistentAssoc03.2.RData"))

#Load in our MBG trait SNPs
load(file.path(data.path, "data_out/BetaCor_05.3.RData"))

#Rename snps file with the new naming
#Specify a function for cleaning directory names
clean_directory_name <- function(x) {
  x %>%
    stringr::str_replace_all("[^[:alnum:]]+", "_") %>% 
    stringr::str_replace_all("_+", "_")   
}


names(snps_df) <- clean_directory_name(names(snps_df))


####################################
####### Load in LD SNPs ############
####################################

#Get the LD SNPs from our covariates of interest - lets start with r >= 0.1
ld_snps_out <- file.path(data.path, "phenosumstats/ldsnps4alphadiv.RData")

covars <- candidate_bugs$prs %>% unique()

if(!file.exists(ld_snps_out)){
  ld_snps <- list()
  for(covar in covars){
    ld_snps[[covar]] <- read.table(file.path(data.path, "phenosumstats", covar, "corsnps", paste0("corsnps_", covar, ".ld")), header=T) %>%
      filter(R2 >= 0.1)
  }
  
  #Save this is an .RData file so we can load locally
  save(ld_snps, file = ld_snps_out)
}else{load(ld_snps_out)}

ld_snps_long <- lapply(names(ld_snps), function(x) {
  cbind(ld_snps[[x]], x) %>% as.data.frame()
}) %>% do.call(rbind, .)


############################################
# Look for enrichment in FGFP summary stats
############################################


#Remove urea as no longer in the final analysis
#ld_snps_long <- filter(ld_snps_long, x != "Serum_urea_levels")

#Focus on SNPs in high LD
ld_snps_long_highr2 <- filter(ld_snps_long, R2 > 0.8)

#Edit naming so nice for display
ld_snps_long_highr2 <- ld_snps_long_highr2 %>%
  mutate(covar_clean = case_when(
    x == "yengo_bmi" ~ "BMI",
    x == "Alanine_aminotransferase_levels" ~ "AAT",
    x == "Gamma_glutamyltransferase_levels_UKB_data_field_30730_" ~ "GGT",
    x == "Hip_circumference_UKB_data_field_49_" ~ "Hip circumference",
    x == "Standing_height" ~ "Standing height",
    x == "Red_blood_cell_count" ~ "RBC count",
    x == "HDL_cholesterol" ~ "HDL",
    x == "Serum_uric_acid_levels" ~ "Uric Acid",
    x == "triglycerides" ~ "Triglycerides",
    x == "Serum_urea_levels" ~ "Urea"
  ))


#Load in data which has all FGFP results with p < 1e-5
load(file.path(data.path, "davesRNTsumstats", "fgfpfilt_sumstats.RData")) 

#These were saved with p < 1e-4, so lets loop through the list and filter out
#Lets collapse into a single df
sigsnps <- setNames(
  lapply(names(out), function(x) {
    filter(out[[x]], frequentist_add_pvalue < 1e-5 & rsid %in% ld_snps_long_highr2$SNP_B)
  }),
  names(out)
)

sigsnps <- sigsnps[sapply(sigsnps, nrow) > 0]

#Print if any SNPs have p < 5e-8
for(i in 1:length(sigsnps)){
  traitname <- names(sigsnps)[i]
  sigsnp_trait <- sigsnps[[i]]
  if(any(sigsnp_trait$frequentist_add_pvalue < 5e-8)){
    cat(paste0("Trait: ", traitname, " has SNPs with p < 5e-8\n"))
  }
}


#Now go through and match the SNPs which were flagged
matched_snps <- list()

#Now we can make a new list which holds the matched microbiome GWAS and covariate SNP data
for(i in 1:length(sigsnps)){
  traitname <- names(sigsnps)[i]
  sigsnp_trait <- sigsnps[[i]]
  
  #Get the covariate SNPs which match
  covar_snp <- filter(ld_snps_long_highr2, SNP_B %in% sigsnp_trait$rsid) %>% pull(SNP_A) %>% unique()
  
  #Make a nested loop incase more than one covar_SNP has been extracted
  for(j in covar_snp){
  
  #Get the covariate name
  covarname <- unique(filter(ld_snps_long_highr2, SNP_A == j)$x)
  
  #Get the covariate SNP data
  covar_snp_data <- snps_df[[covarname]] %>%
    filter(rsid %in% j)
  
  #Get the microbiome GWAS SNP data
  fgfp_snp_data <- sigsnp_trait %>%
    filter(rsid %in% filter(ld_snps_long_highr2, SNP_A == j)$SNP_B)
  
  #Store in list
  matched_snps[[paste0(traitname, "_", j)]] <- list(
    "covar_snp_data" = covar_snp_data,
    "fgfp_snp_data" = fgfp_snp_data,
    "covarname" = covarname,
    "R2" = unique(filter(ld_snps_long_highr2, SNP_A == j & SNP_B %in% fgfp_snp_data$rsid)))
  }
}  

  # Initialize results list
  summary_table <- list()
  
  for (k in seq_along(matched_snps)) {
    
    snpofinterest <- matched_snps[[k]]$fgfp_snp_data %>%
      mutate(
        # Identify palindromic SNPs
        is_palindromic = (alleleB == "A" & alleleA == "T") |
          (alleleB == "T" & alleleA == "A") |
          (alleleB == "C" & alleleA == "G") |
          (alleleB == "G" & alleleA == "C")) %>%
      filter(frequentist_add_pvalue == min(frequentist_add_pvalue))
    
      
    
    # --- Combine all info into one tibble row ---
    summary_table[[k]] <- tibble(
      Micro_Trait = gsub("(_RNT|_HB).*", "", names(matched_snps)[k]),
      Covar = matched_snps[[k]]$covarname,
      covar_clean = unique(matched_snps[[k]]$R2$covar_clean),
      Micro_rsid = snpofinterest$rsid[1],
      Micro_snpid = snpofinterest$snpid[1],
      Micro_chr = snpofinterest$chromosome[1],
      Micro_bp = snpofinterest$position[1],
      Micro_beta = snpofinterest$frequentist_add_beta_1[1],
      Micro_se = snpofinterest$frequentist_add_se_1[1],
      Micro_p = snpofinterest$frequentist_add_pvalue[1],
      Micro_RA = snpofinterest$alleleA[1],
      Micro_EA = snpofinterest$alleleB[1],
      Micro_EAF = snpofinterest$all_maf[1],
      is_palindromic = snpofinterest$is_palindromic[1],
      
      Covar_rsid = matched_snps[[k]]$covar_snp_data$rsid[1],
      Covar_snpid = matched_snps[[k]]$covar_snp_data$fgfp_linker[1],
      Covar_beta = matched_snps[[k]]$covar_snp_data$beta[1],
      Covar_se = matched_snps[[k]]$covar_snp_data$SE[1],
      Covar_p = matched_snps[[k]]$covar_snp_data$pval[1],
      Covar_EA = matched_snps[[k]]$covar_snp_data$eff.allele[1],
      Covar_EAF = matched_snps[[k]]$covar_snp_data$EAF[1],
      R2 =  filter(matched_snps[[k]]$R2, SNP_B == snpofinterest$rsid)[["R2"]]
    ) %>% as.data.frame()
  }
  
# --- Combine all rows into a single data frame ---
summary_table_fgfp <- bind_rows(summary_table)



############################################
# Now do MBG
############################################

  
  
#Need to decide how to harmonise SNPs off the covariate one or the Microbiome SNP?  

#Load in the MGBe4_snps
load(file.path(data.path, "mibiogen_ss/MBGe4_snps.RData"))
mbg_snps <- filter(MBGe4_snps, P.weightedSumZ < 1e-5)

#See if any of the covariate SNPs of interest are here
sigsnps_mbg <- mbg_snps %>%
  filter(rsID %in% ld_snps_long_highr2$SNP_B)

#Are any of these genome-wide significant?
#sigsnps_mbg <- filter(sigsnps_mbg, P.weightedSumZ < 5e-8)

#Make into a list split by "bac" column so code will run
sigsnps_mbg <- split(sigsnps_mbg, sigsnps_mbg$bac)

#Again with how SNPs were matched to covariate with FGFP, do the same with MiBioGen
matched_snps_mbg <- list()

for(i in 1:length(sigsnps_mbg)){
  traitname <- names(sigsnps_mbg)[i]
  sigsnp_trait <- sigsnps_mbg[[i]]
  
  #Get the covariate SNPs which match
  covar_snp <- filter(ld_snps_long_highr2, SNP_B %in% sigsnp_trait$rsID) %>% pull(SNP_A) %>% unique()
  
  #Make a nested loop incase more than one covar_SNP has been extracted
  for(j in covar_snp){
    
    #Get the covariate name
    covarname <- unique(filter(ld_snps_long_highr2, SNP_A == j)$x)
    
    #Get the covariate SNP data
    covar_snp_data <- snps_df[[covarname]] %>%
      filter(rsid %in% j)
    
    #Get the microbiome GWAS SNP data
    mbg_snp_data <- sigsnp_trait %>%
      filter(rsID %in% filter(ld_snps_long_highr2, SNP_A == j)$SNP_B)
    
    #Store in list
    matched_snps_mbg[[paste0(traitname, "_", j)]] <- list(
      "covar_snp_data" = covar_snp_data,
      "mbg_snp_data" = mbg_snp_data,
      "covarname" = covarname,
      "R2" = unique(filter(ld_snps_long_highr2, SNP_A == j & SNP_B %in% mbg_snp_data$rsID)))
  }
}
  
#Now make a clean df of this

# Initialize results list
summary_table_mbg <- list()
for (k in seq_along(matched_snps_mbg)) {
  
  snpofinterest <- matched_snps_mbg[[k]]$mbg_snp_data %>%
    mutate(
      # Identify palindromic SNPs
      is_palindromic = (eff.allele == "A" & ref.allele == "T") |
        (eff.allele == "T" & ref.allele == "A") |
        (eff.allele == "C" & ref.allele == "G") |
        (eff.allele == "G" & ref.allele == "C")) %>%
    filter(P.weightedSumZ == min(P.weightedSumZ))
  
  
  
  # --- Combine all info into one tibble row ---
  summary_table_mbg[[k]] <- tibble(
    Micro_Trait = gsub("\\.id\\..*", "", names(matched_snps_mbg)[k]),
    Covar = matched_snps_mbg[[k]]$covarname,
    covar_clean = unique(matched_snps_mbg[[k]]$R2$covar_clean),
    Micro_rsid = snpofinterest$rsID[1],
    Micro_chr = snpofinterest$chr[1],
    Micro_bp = snpofinterest$bp[1],
    Micro_beta = snpofinterest$beta[1],
    Micro_se = snpofinterest$SE[1],
    Micro_p = snpofinterest$P.weightedSumZ[1],
    Micro_EA = snpofinterest$eff.allele[1],
    Micro_RA = snpofinterest$ref.allele[1],
    #Micro_EAF = snpofinterest$EAF[1],
    is_palindromic = snpofinterest$is_palindromic[1],
    
    Covar_rsid = matched_snps_mbg[[k]]$covar_snp_data$rsid[1],
    Covar_beta = matched_snps_mbg[[k]]$covar_snp_data$beta[1],
    Covar_se = matched_snps_mbg[[k]]$covar_snp_data$SE[1],
    Covar_p = matched_snps_mbg[[k]]$covar_snp_data$pval[1],
    Covar_EA = matched_snps_mbg[[k]]$covar_snp_data$eff.allele[1],
    Covar_EAF = matched_snps_mbg[[k]]$covar_snp_data$EAF[1],
    R2 =  filter(matched_snps_mbg[[k]]$R2, SNP_B == snpofinterest$rsID)[["R2"]]
  ) %>% as.data.frame()
}

#Make to a table
summary_table_mbg <- bind_rows(summary_table_mbg)
    
# Now can output these results to make some locus plots in a new script
save(summary_table_fgfp, summary_table_mbg, file = file.path(data.path, "data_out/Enrichment_SNPs_FGFP_MBG.RData"))







