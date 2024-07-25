####################################################################
#This script will run our regressions of interest
####################################################################

#Set wd to 03_Regressions !

####################################################################
#First run the covariate PRS -> PRS and save lm output
####################################################################

#Set up
#Library
library(dplyr)
library(data.table)
library(broom)

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

#Extract the PRSs we want to use, removed duplicates and any factor variables at the moment
prs2use <- data4reg[data4reg$pheno_vartype == "cont", "prs"] %>% unique()

#Lets make df for the scaled data so can also have these coefs if needed
prspheno_df_scale <- scale(data4prs$pheno_covariate_prs[,c(prs2use, contvars, confounders)]) %>% as.data.frame()

#Now lets extract our coefficents from the model
###

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
m_scale <- lm(reformulate(c(confounders, pheno), response = prs), data = prspheno_df_scale)
out2 <- tidy(m_scale, conf.int = T) %>% dplyr::filter(term == pheno) %>% cbind(nobs(m_scale))
regout_covarprs_cont_scale <- rbind(regout_covarprs_cont_scale, out2) 
}


#Combine the PRS name with output
regout_covarprs_cont <- cbind(regout_covarprs_cont, data4reg[data4reg$pheno_vartype == "cont", "prs"])
regout_covarprs_cont_scale <- cbind(regout_covarprs_cont_scale, data4reg[data4reg$pheno_vartype == "cont", "prs"])



#Save the output to object
regout_covarprs <- list(data4reg = data4reg, 
     regout_covarprs_cont = regout_covarprs_cont,
     regout_covarprs_cont_scale = regout_covarprs_cont_scale)



####################################################################
#Next we can run the linear regression of PRS -> RNT MT
####################################################################

#First lets make a df with all our variables of interest, this will be our confounders, cont PRSs and MTs
prsmicropheno_df <- merge(data4prs$pheno_covariate_prs[,c(prs2use, contvars, "fgfp_id")],
                     dplyr::select(data4prs$gwasedmts, linker, ends_with("RNTRes")),
                     by.x = "fgfp_id",
                     by.y = "linker",
                     all = T) %>% distinct()

#Also make a scaled df:
prsmicropheno_df_scale <- scale(prsmicropheno_df[,-1]) %>% as.data.frame()

#Extract the RNT traits
RNTs <- data4prs$gwasedmts %>% dplyr::select(ends_with("RNTRes")) %>% colnames()

#Empty tibbles for output of results
#We can run this both as a multivariate and univariate regression to see how results compare
regout_RNTs_univar <- tibble()
regout_RNTs_multivar <- tibble()

for(mt in RNTs){
  #First run multivar
  m <- lm(reformulate(prs2use, response = mt), data = prsmicropheno_df_scale)
  out <- tidy(m, conf.int = TRUE) %>% cbind(nobs(m), mt) %>% dplyr::filter(term != "(Intercept)")
  regout_RNTs_multivar <- rbind(regout_RNTs_multivar, out) 
  
  #Then run it as univar
  for (prs in prs2use){
 
#Run the regression, scale it, extract coefficents of interest and bind to df
    lm2 <- lm(reformulate(prs, response = mt), data = prsmicropheno_df_scale)
    out2 <- tidy(lm2) %>% 
      dplyr::filter(term == prs) %>% 
      cbind(nobs(lm2), mt)
  
  #Save output
  regout_RNTs_univar <- rbind(regout_RNTs_univar, out2) 
  }
}
  
#We can not save our output
regout_prsRNTs <- list(regout_RNTs_univar = regout_RNTs_univar,
                       regout_RNTs_multivar = regout_RNTs_multivar)



####################################################################
#Next we can run the observational regression of covar -> RNT MT
####################################################################
#Load in h2 bugs object
load(file = file.path(data.path, "data_out/h2bugs_02.2.RData"))

#For these traits with PRS signal do we also see an observational signal?
signal <- lapply(regout_prsRNTs, 
                 function(x) filter(x, mt %in% paste0(RNT_h2_traits$TaxaName, "_RNTRes") &
                                      p.value < 0.05))

#Can just look at those with signal as highlighted above, however I think best to run for all then filter on the ones we are interested in after

obs_est <- tibble()
for(mt in RNTs){

for(pheno in contvars)  {
    lm <- lm(reformulate(pheno, response = mt), data = prsmicropheno_df_scale)
    out <- lm %>% 
      tidy() %>% 
      dplyr::filter(term != "(Intercept)") %>% 
      cbind(nobs(lm), mt)
    
    obs_est <- rbind(obs_est, out) 
}
}


##########################################################################################
#Now we can save all objects to one .RData object so easier to load into the Markdown
##########################################################################################
save(regout_covarprs, regout_prsRNTs, obs_est,
     file = file.path(data.path, "data_out/data_Reg03.1.RData"))

