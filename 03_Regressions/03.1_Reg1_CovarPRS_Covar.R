####################################################################
#This step will regress the covariate PRS -> PRS and save lm output
####################################################################

#Set wd to 03_Regressions !

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
load(file.path(data.path, "data_out/fgfpdata4prs.RData"))

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
              "e.GFR")

#Break the df into continous and factor label
data4reg$pheno_vartype <- ifelse(data4reg$pheno %in% contvars, "cont", "fact")

#Now lets extract our coefficents from the model
#I have not standardised these variables and so are in the scale of the original GWAS

regout_covarprs_cont <- tibble()
regout_covarprs_cont_scale <- tibble()
for (i in 1:sum(data4reg$pheno_vartype == "cont")){

#Extract name of pheno and prs
pheno <- data4reg[data4reg$pheno_vartype == "cont", "pheno"][i]
prs <- data4reg[data4reg$pheno_vartype == "cont", "prs"][i]

#First run regression unscaled
m <- lm(reformulate(c(confounders, pheno), response = prs), data = data4prs$pheno_covariate_prs)
out <- tidy(m, conf.int = TRUE) %>% dplyr::filter(term == pheno) %>% cbind(nobs(m))
regout_covarprs_cont <- rbind(regout_covarprs_cont, out) 

#Then scale
m_scale <- lm.beta(m)
out2 <- tidy(m_scale) %>% dplyr::filter(term == pheno) %>% cbind(nobs(m))
regout_covarprs_cont_scale <- rbind(regout_covarprs_cont_scale, out2) 
}

#Fix the CIs on the scaled
regout_covarprs_cont_scale <- regout_covarprs_cont_scale %>%
  mutate(conf.low = std_estimate - 1.96 * std.error,
         conf.high = std_estimate + 1.96 * std.error)

#Combine the PRS name with output
regout_covarprs_cont <- cbind(regout_covarprs_cont, data4reg[data4reg$pheno_vartype == "cont", "prs"])
regout_covarprs_cont_scale <- cbind(regout_covarprs_cont_scale, data4reg[data4reg$pheno_vartype == "cont", "prs"])



#Save the output
#This will be a list with the data4reg file and regression coefficents
#Still need to add the factor vars and the standardised regressions if going forward with these analyses
regout_covarprs <- list(data4reg = data4reg, 
     regout_covarprs_cont = regout_covarprs_cont,
     regout_covarprs_cont_scale = regout_covarprs_cont_scale)

save(regout_covarprs, 
     file = file.path(data.path, "data_out/regout_covarprs.RData"))



