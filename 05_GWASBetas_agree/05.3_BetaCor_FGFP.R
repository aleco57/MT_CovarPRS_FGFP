#### This script will look at our candidate relationships and see if these covariate SNPs are showing in the sumstats for those specific traits in mibiogen

#Set up
#Library
library(tidyverse)
library(data.table)
library(ggrepel)
library(broom)
library(TwoSampleMR)


#Data path for proj dir
data.path <- "../../../data/PhenoPRS"

#Also load in path to scratch, proj_dir and fgfp data
source("../../parameters/base_dir.R")
source("../../parameters/scratch_path.txt")

#Load in the bugs we are interested in 
load(file.path(data.path, "data_out/data_ConsistentAssoc03.2.RData"))

#Also load in the SNP data
load(file.path(data.path, "phenosumstats/filtered_sumstats_list.RData"))

load(file.path(data.path, "data_out/filtered_bug_sumstats_05.1.RData"))


#Lets look at FGFP sumstats
############################


#Get our candidiates we are intersted in
#

harm_data_all <- list()
lm_out <- c()
for(covar in names(merged_gwasbetas)){
  for(bug in names(merged_gwasbetas[[covar]])){
    
    #Harmonise so same effect allele the compute the correlation between the betas
    fgfp_bug <- merged_gwasbetas[[covar]][[bug]] %>%
      data.frame() %>%
      format_data(type = "exposure",
                  snp_col = "rsid.x",
                  beta_col = "frequentist_add_beta_1",
                  se_col = "frequentist_add_se_1",
                  eaf_col = "all_maf",
                  effect_allele_col = "alleleB",
                  other_allele_col = "alleleA")
    
    if(covar == "yengo_bmi"){
      covarss <- merged_gwasbetas[[covar]][[bug]][,26:32] %>%
        data.frame() %>%
        format_data(type = "outcome",
                    snp_col = "rsid.y",
                    beta_col = "beta",
                    se_col = "SE",
                    eaf_col = "EAF",
                    effect_allele_col = "eff.allele",
                    other_allele_col = "Oth_al")
    }else {
      covarss <- merged_gwasbetas[[covar]][[bug]][,26:31] %>%
        data.frame() %>%
        format_data(type = "outcome",
                    snp_col = "rsid.y",
                    beta_col = "beta",
                    se_col = "SE",
                    eaf_col = "EAF",
                    effect_allele_col = "eff.allele")
    }
    
    
    harm_data_all[[covar]][[bug]] <- harmonise_data(exposure = fgfp_bug, outcome = covarss)
    
    #Linear regression
    out <- lm(beta.exposure ~ beta.outcome, data = harm_data_all[[covar]][[bug]]) %>% summary() %>% tidy() %>% filter(term == "beta.outcome")
    
    #We can also run IVW regression so we can uncertainty of our SNP estimates into the model
    harm_data_all[[covar]][[bug]]$weights <- 1/(harm_data_all[[covar]][[bug]]$se.exposure^2 + harm_data_all[[covar]][[bug]]$se.outcome^2)
    ivw_out <- lm(beta.exposure ~ beta.outcome, data = harm_data_all[[covar]][[bug]], weights = harm_data_all[[covar]][[bug]]$weights) %>% summary() %>% tidy() %>% filter(term == "beta.outcome")
    colnames(ivw_out) <- paste0("ivw_", colnames(ivw_out))
    
    #Now we can merge with our estimates
    out <- cbind(out, ivw_out, covar, bug)
    lm_out <- rbind(lm_out, out)
  }}

#Save the output
save(lm_out, harm_data_all, file = file.path(data.path, "data_out/BetaCor_05.3.RData"))


