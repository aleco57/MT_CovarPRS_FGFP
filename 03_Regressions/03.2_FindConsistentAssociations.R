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

data4prs$matchedvars$prslabel <- paste0(data4prs$matchedvars$prs, ".Pt_5e.08")

#Merge the phenotype information so we can match prs and observational results
prs2merge <- merge(regout_prsRNTs$regout_RNTs_univar,
                   data4prs$matchedvars,
                   by.x = "term",
                   by.y="prslabel")

#Merge PRS and obs results, only those with obs result with p < 0.05
obsest_prsest <- merge(obs_est,
                       prs2merge,
                       by.x = c("term", "mt"),
                       by.y = c("pheno", "mt"),
                       all.x = T)

#Edit colnames so we know which beta belongs to which regression
colnames(obsest_prsest) <- gsub(".x", "_obs", colnames(obsest_prsest))
colnames(obsest_prsest) <- gsub(".y", "_prs", colnames(obsest_prsest))

#Where do we have best candidates for microbial trait just capturing the covariate?
candidate_bugs <- filter(obsest_prsest, p.value_obs < 0.05 & p.value_prs < 0.05 & estimate_obs * estimate_prs > 0)

#Save the output so can be used for subsequent analyses
save(obsest_prsest, candidate_bugs,
     file = file.path(data.path, "data_out/data_ConsistentAssoc03.2.RData"))


