#######################################################################
##### Which Microbial traits show p < 0.05 observational association 
##### AND directionality agreement for PRS analysis?
#######################################################################

###########
# Set up
###########
library(dplyr)
library(readxl)

#Data path for proj dir
data.path <- "../../../data/PhenoPRS"

#Load in our data
load(file.path(data.path, "data_out/fgfpdata4prs.RData"))

#Read in the regressions output so we can see which bug is associated with which PRS
load(file.path(data.path, "data_out/data_Reg03.1.RData"))

#Removing unwanted PRSs
data4reg <- data4prs[["matchedvars"]] %>% 
  dplyr::filter(!prs %in% c("stool_freq.Pt_5e.08", "creatinine3.Pt_5e.08", "creatinine2.Pt_5e.08"))

#Merge the phenotype information so we can match prs and observational results
prs2merge <- merge(regout_prsRNTs$regout_RNTs_univar,
                   data4prs$matchedvars,
                   by.x = "term",
                   by.y="prs")

#Merge PRS and obs results, only those with obs result with p < 0.05
obsest_prsest <- merge(obs_est[obs_est$p.value<0.05,],
                       prs2merge,
                       by.x = c("term", "mt"),
                       by.y = c("pheno", "mt"),
                       all.x = T) %>%
  #And remove creatinine for now as may have to generate a new PRS for this
  filter(!term %in% c("Creatinine_mgdL", "e.GFR"))

#Filter to those where we have beta agreement
obsest_prsest <- obsest_prsest[obsest_prsest$estimate.x * obsest_prsest$estimate.y > 0,]


#Edit colnames so we know which beta belongs to which regression
colnames(obsest_prsest) <- gsub(".x", "_obs", colnames(obsest_prsest))
colnames(obsest_prsest) <- gsub(".y", "_prs", colnames(obsest_prsest))

#Save the output so can be used for subsequent analyses
save(obsest_prsest,
     file = file.path(data.path, "data_out/data_ConsistentAssoc03.2.RData"))
