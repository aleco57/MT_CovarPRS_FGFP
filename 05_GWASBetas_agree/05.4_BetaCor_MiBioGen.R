#### Beta agree for MiBioGen GWAS results
### What would be more informative is if we observe beta coefficients agreeing in the independent MiBioGen dataset

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

load(file.path(data.path, "data_out/filtered_bug_sumstatsMBG_05.1.RData"))


#Make matched mbg_fgfp df
## First match our traits of interest with the mibiogen SNPs

#Make our matched mbg_fgfp bug dataframe
mbgbugs <- c("simpson", #Div_NumGen
             "invsimpson", #Div_Chao
             "family.Lachnospiraceae.id.1987",
             "genus.Barnesiella.id.944",
             "genus.Roseburia.id.2012", 
             "family.Ruminococcaceae.id.2050", #G_Sporobacter
             "family.Ruminococcaceae.id.2050", #G_unclassified
             NA, #G_Unclassified__K_Bacteria
             "order.Clostridiales.id.1863",
             "phylum.Firmicutes.id.1672", #G_Unclassified
             "family.Erysipelotrichaceae.id.2149",
             "genus.Holdemanella.id.11393",
             "genus.Ruminococcus1.id.11373",
             "family.Erysipelotrichaceae.id.2149", #G_unclassified
             "family.Ruminococcaceae.id.2050",
             "genus.Oscillibacter.id.2063",
             "genus.Bacteroides.id.918",
             "family.Rhodospirillaceae.id.2717", #G_Aestuariispira
             "family.Lachnospiraceae.id.1987", #G_Coprococcus
             "shannon",
             "family.Porphyromonadaceae.id.943")

matched_bugs <- data.frame(fgfpbugs = unique(candidate_bugs$mt), mbgbugs = mbgbugs)

#Remove the NA
matched_bugs <- matched_bugs[!is.na(matched_bugs$mbgbugs),]

harm_data_all <- list()
lm_out <- c()
for(covar in names(merged_gwasbetas_mbg)){
  for(bug in names(merged_gwasbetas_mbg[[covar]])){
    
    #Harmonise so same effect allele the compute the correlation between the betas
    fgfp_bug <- merged_gwasbetas_mbg[[covar]][[bug]] %>%
      data.frame() %>%
      format_data(type = "exposure",
                  snp_col = "rsID",
                  beta_col = "beta.x",
                  se_col = "SD",
                  eaf_col = "weighted_MAF",
                  effect_allele_col = "eff.allele.x",
                  other_allele_col = "ref.allele",
                  pval_col = "P.weightedSumZ")
    
    if(covar == "yengo_bmi"){
      covarss <- merged_gwasbetas_mbg[[covar]][[bug]] %>%
        data.frame() %>%
        format_data(type = "outcome",
                    snp_col = "rsID",
                    beta_col = "beta.y",
                    se_col = "SE",
                    eaf_col = "EAF",
                    effect_allele_col = "eff.allele.y",
                    other_allele_col = "Oth_al",
                    pval_col = "pval")
    }else {
      covarss <- merged_gwasbetas_mbg[[covar]][[bug]] %>%
        data.frame() %>%
        format_data(type = "outcome",
                    snp_col = "rsID",
                    beta_col = "beta.y",
                    se_col = "SE",
                    eaf_col = "EAF",
                    effect_allele_col = "eff.allele.y",
                    pval_col = "pval")
    }
    
    
    harm_data_all[[covar]][[bug]] <- harmonise_data(exposure = fgfp_bug, outcome = covarss)
    
    #Linear regression
    out <- lm(beta.exposure ~ beta.outcome, data = harm_data_all[[covar]][[bug]]) %>% summary() %>% tidy() %>% filter(term == "beta.outcome")
    
    #IVW Linear Regression
    #We can weight by the exposure and outcome
    harm_data_all[[covar]][[bug]]$weights <- 1/(harm_data_all[[covar]][[bug]]$se.exposure^2 + harm_data_all[[covar]][[bug]]$se.outcome^2)
    
    ivw_out <- lm(beta.exposure ~ beta.outcome, data = harm_data_all[[covar]][[bug]], weights = harm_data_all[[covar]][[bug]]$weights) %>% summary() %>% tidy() %>% filter(term == "beta.outcome")
    colnames(ivw_out) <- paste0("ivw_", colnames(ivw_out))
    out <- cbind(out, ivw_out, covar, bug)
    lm_out <- rbind(lm_out, out)
  }}


#Save the output
lm_out_mbg <- lm_out %>% rename(mbgbugs = bug)
harm_data_all_mbg <- harm_data_all

# Add in the matched bugs to lm
lm_out_mbg <- right_join(lm_out_mbg,
                         matched_bugs,
                         by = "mbgbugs")

#Filter out the wrong matched

#First split the data into covars so we can filter out the incorrect pairings
lm_out_mbg_split <- split(lm_out_mbg, lm_out_mbg$covar)

#Now for each covar check that the macthed bug is correct
for(covar in names(lm_out_mbg_split)){
  
  #Get the right bugs
  bugs <- filter(candidate_bugs, prs == covar)[["mt"]]
  
  lm_out_mbg_split[[covar]] <-  filter(lm_out_mbg_split[[covar]], fgfpbugs %in% bugs)
}


lm_out_mbg <- do.call(rbind, lapply(names(lm_out_mbg_split), function(covar) {
  data <- lm_out_mbg_split[[covar]]
  data$covar <- covar  # Add the list element name as a column
  data
}))

rownames(lm_out_mbg) <- NULL  # Remove rownames

save(lm_out_mbg, harm_data_all_mbg, file = file.path(data.path, "data_out/BetaCorMBG_05.3.RData"))






