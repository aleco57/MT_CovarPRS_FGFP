####################################################################
#This step will regress the covariate PRS -> PRS and save lm output
####################################################################

#Set up
#Library
library(dplyr)
library(data.table)
library(lm.beta)

#Data path for proj dir
data.path <- "../../../data/PhenoPRS"

#Also load in path to scratch, proj_dir and fgfp data
source("../../parameters/base_dir.R")
gwas_table <- read.csv(file = file.path(data.path, "ProcessedGWASTable_OG_GC.csv"), header = T)
load(file.path(data.path, "data4h2rmd/fgfpdata4prs.RData"))

#Also make a vector for our confounders of interest
PCs <- c()
for (i in 1:10) {
  PCs <- c(PCs, paste0("PC", i))
}
confounders <- c(PCs, "age", "snpsex")
confoundbatch <- c(PCs, "age", "snpsex", "drilled.cut", "Alliquoted.by", "Aliquote.date_year", "Date.of.extraction_year", "PCR.Plate")

### Now edit our "data4prs[["matchedvars"]]" slightly to remove PRSs we are no longer interested in
data4reg <- data4prs[["matchedvars"]] %>% 
  dplyr::filter(!prs %in% c("stool_freq.Pt_5e.08", "creatinine3.Pt_5e.08"))

##Regress our covar PRS -> covar
#I don't think you should be running a regression model with the factor variables

#We are going to start our pipeline with just the continuous variables
#Have also left sleeping hours per day in as it is nicely normally distributed
#Only vars which have been left out are the "average consumption" vars
contvars <- c("Hemoglobine_gdL",
              "RBC_milj_mm3",
              "Urinezuur_mgdL",
              "HDL.chol_mgdL",
              "BMI",
              "Creatinine_mgdL",
              "Triglyceriden_mgdL",
              "e.GFR",
              "sleeping_hours_per_day")

#Now lets extract our coefficents from the model
#I have not standardised these variables and so are in the scale of the original GWAS
regout_covarprs_cont <- tibble()
for (i in 1:nrow(data4reg)){
m <- lm(reformulate(c(confounders, data4reg$pheno[i]), response = data4reg$prs[i]), data = data4prs$pheno_covariate_prs)
#m_scale <- lm.beta(m)
out <- tidy(m, conf.int = TRUE) %>% dplyr::filter(term == data4reg$pheno[i]) %>% cbind(nobs(m))
regout_covarprs_cont <- rbind(regout_covarprs_cont, out) 
}

#Combine the PRS name with output
regout_covarprs <- cbind(regout_covarprs, data4reg$prs)

#Save the output
save(regout_covarprs, file = file.path(data.path, "data4h2rmd/regout_covarprs.RData"))



