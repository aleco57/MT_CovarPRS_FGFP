#### For our candidate bugs - lets also look if this relationship is repliacted in LogRes transformation so we can take these forwards to MiBioGen sumstats

#Set up
#Library
library(tidyverse)
library(data.table)
library(broom)
library(sensemakr)
library(ggplot2)


#Data path for proj dir
data.path <- "../../../data/PhenoPRS"

#Also load in path to scratch, proj_dir and fgfp data
source("../../parameters/base_dir.R")
load(file.path(data.path, "data_out/fgfpdata4prs.RData"))

#Also make a vector for our confounders of interest
PCs <- c()
for (i in 1:10) {
  PCs <- c(PCs, paste0("PC", i))
}
confounders <- c(PCs, "age", "snpsex")
confoundbatch <- c(PCs, "age", "snpsex", "drilled.cut", "Alliquoted.by", "Aliquote.date_year", "Date.of.extraction_year", "PCR.Plate")


#Now load the data
load(file.path(data.path, "data_out/data_ConsistentAssoc03.2.RData"))



##### Now run observational and PRS on LogRes traits
logmt_names <- gsub("_RNTRes", "_LOGRes", candidate_bugs$mt)
data4prs[["gwasedmts"]] %>% select()

#Make df
prs2use <- unique(candidate_bugs$term_prs)

#First lets make a df with all our variables of interest, this will be our confounders, cont PRSs and MTs
prsmicropheno_df <- merge(data4prs$pheno_covariate_prs[,c(data4prs$phenos2use$pheno, prs2use, "fgfp_id")],
                          dplyr::select(data4prs$gwasedmts, linker, all_of(logmt_names)),
                          by.x = "fgfp_id",
                          by.y = "linker",
                          all = T) %>% distinct()
#Also make a scaled df:
prsmicropheno_df_scale <- scale(prsmicropheno_df[,-1]) %>% as.data.frame()

#Empty tibbles for output of results
#We can run this both as a multivariate and univariate regression to see how results compare
regout_LOGs_univar <- tibble()

for(i in 1:nrow(candidate_bugs)){
    mt <- logmt_names[i]
    prs <- candidate_bugs$term_prs[i]
    #Run the regression, scale it, extract coefficents of interest and bind to df
    lm2 <- lm(reformulate(prs, response = mt), data = prsmicropheno_df_scale)
    out2 <- tidy(lm2) %>% 
      dplyr::filter(term == prs) %>% 
      cbind(nobs(lm2), mt)
    
    #Save output
    regout_LOGs_univar <- rbind(regout_LOGs_univar, out2) 
}

regout_LOGs_univar$RNTlink <- candidate_bugs$mt
colnames(regout_LOGs_univar)[2:5] <- paste0(colnames(regout_LOGs_univar)[2:5], "LogRes")

# Lets merge with the candidate bugs df
regout_LOGs_univar <- merge(regout_LOGs_univar, candidate_bugs, by.x = c("term", "RNTlink"), by.y = c("term_prs", "mt"))

#Do signs all agree?
sign(regout_LOGs_univar$estimateLogRes) == sign(regout_LOGs_univar$estimate_prs)

#Yes they do ! Happy to take forward to look at mibiogen sumstats

#Saveoutput
save(regout_LOGs_univar,
     file = file.path(data.path, "data_out/data_RegLog03.3.RData"))
