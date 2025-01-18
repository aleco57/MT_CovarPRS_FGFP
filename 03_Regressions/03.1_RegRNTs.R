####################################################################
#This script will run our regressions of interest
####################################################################

#Set wd to 03_Regressions !

####################################################################
#First run the covariate PRS -> PRS and save lm output
####################################################################

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
age_sex <- c("age", "snpsex")
### Now edit our "data4prs[["matchedvars"]]" slightly so we can select the PRSs we would like to take forwards
#First we want to label if our variables are cont / binary / ordinal
#We will use the following labels:
# 1 = binary
# 2 = ordinal
# 3 = continuous

#See how many unique obs there are, if 3 we will label as binary, if more than 10 then we will label as continuous
num_obs <- apply(select(data4prs$pheno_covariate_prs, data4prs$matchedvars$pheno), 2, function(x) length(unique(x)))

data4prs$matchedvars$vartype <- ifelse(data4prs$matchedvars$pheno %in% names(num_obs[num_obs == 3]), 1, 2)
data4prs$matchedvars$vartype <- ifelse(data4prs$matchedvars$pheno %in% names(num_obs[num_obs > 9]), 3, data4prs$matchedvars$vartype)

##Regress our covar PRS -> covar
#Lets start with the continuous variables

#Lets make df for the scaled data so can also have these coefs if needed
#prspheno_df_scale <- scale(data4prs$pheno_covariate_prs[,c(prs2use, contvars, confounders)]) %>% as.data.frame()

#Now lets extract our coefficents from the model and R2

#regout_covarprs_cont_scale <- tibble()


#Run the loop for continuous variables
regout_covarprs_cont <- tibble()
for (i in which(data4prs$matchedvars$vartype == 3)){

#Extract name of pheno and prs
pheno <- data4prs$matchedvars[i, "pheno"]
prs <- paste0(data4prs$matchedvars[i, "prs"], ".Pt_5e.08")

#First run regression unscaled
m <- lm(reformulate(c(confounders, prs), response = pheno), data = data4prs$pheno_covariate_prs)
out <- tidy(m, conf.int = TRUE) %>% dplyr::filter(term == prs) %>% cbind(nobs(m))

#We also want to add the partial R2 of the model
out <- cbind(out, partial_r2(m)[prs])

regout_covarprs_cont <- rbind(regout_covarprs_cont, out) 
}


##### Which beta coefficients are not significant ?
bad_prs_cont <- gsub(".Pt_5e.08", "", filter(regout_covarprs_cont, p.value > 0.01)[["term"]])


#Now lets look at the binary variables
#########################################

#Now we can look at our binary variables - these will need to be cleaned first
binary_variables <- data4prs$matchedvars[data4prs$matchedvars$vartype == 1, "pheno"]

#Clean binary variables so are all 1/0
data4prs$pheno_covariate_prs[binary_variables] <- 
  data4prs$pheno_covariate_prs[binary_variables] %>%
  mutate(
    # Convert TRUE/FALSE to 0/1
    across(where(is.logical), ~ as.numeric(.)),
    
    # Convert yes/no to 1/0
    across(where(~ is.character(.) && all(. %in% c("yes", "no", NA))), 
           ~ as.numeric(factor(., levels = c("no", "yes"))) - 1),
    
    # If already 0/1, leave as is (no conversion needed for numeric 0/1)
    across(where(is.numeric), ~ .)
  )


#Now we can run logistic regression

regout_covarprs_bin <- tibble()
for (i in which(data4prs$matchedvars$vartype == 1)){
  
  #Extract name of pheno and prs
  pheno <- data4prs$matchedvars[i, "pheno"]
  prs <- paste0(data4prs$matchedvars[i, "prs"], ".Pt_5e.08")
  
  #First run regression unscaled
  m <- glm(reformulate(c(confounders, prs), response = pheno), data = data4prs$pheno_covariate_prs, family = binomial(link = "logit"))
  out <- tidy(m, conf.int = TRUE) %>% dplyr::filter(term == prs) %>% cbind(nobs(m))
  
  regout_covarprs_bin <- rbind(regout_covarprs_bin, out) 
}


#### None of the logistic regression models suggest evidence of association with PRS


#### All of the others are also unlikely to associate, so we will stick with the continuous variables !

phenos2use <- filter(data4prs$matchedvars, vartype == 3 & !prs %in% bad_prs_cont)

data4prs[["phenos2use"]] <- filter(data4prs$matchedvars, vartype == 3 & !prs %in% bad_prs_cont)

data4prs[["r2ofcovars"]] <- regout_covarprs_cont

#Combine the PRS name with output
#regout_covarprs_cont <- cbind(regout_covarprs_cont, data4reg[data4reg$pheno_vartype == "cont", "prs"])
#regout_covarprs_cont_scale <- cbind(regout_covarprs_cont_scale, data4reg[data4reg$pheno_vartype == "cont", "prs"])

#Now we can save our new data object
save(data4prs, file = file.path(data.path, "data_out/fgfpdata4prs.RData"))


####################################################################
#Next we can run the linear regression of PRS -> RNT MT
####################################################################

prs2use <- paste0(data4prs$phenos2use$prs, ".Pt_5e.08")

#First lets make a df with all our variables of interest, this will be our confounders, cont PRSs and MTs
prsmicropheno_df <- merge(data4prs$pheno_covariate_prs[,c("fgfp_id", "IID", data4prs$phenos2use$pheno, prs2use, age_sex)],
                     dplyr::select(data4prs$gwasedmts, linker, ends_with("RNTRes")),
                     by.x = "fgfp_id",
                     by.y = "linker",
                     all = T) %>% distinct()



#We can make a prsmicro scaled removing the genetic IDs we do not want
#First load in the exclusion file - FScalePC does not need to be removed - it is explained in Dave's QC file
excl <- read.table(file.path(pdir, "data/FGFP_data/FGFP_MT/sample_qc/Sample_Exclusion_Criteria.txt"),
                   header = TRUE, as.is = TRUE, sep = "\t")

#This function removes the genetic IDs we do not want and also removes those without microbiome data
IDs2remove <- apply(excl[, c(3, 5:11)], 2, function(column) {
  rownames(excl)[which(column == 1)]
}) %>% unlist() %>% unname() %>% unique()


prsmicro_df <- prsmicropheno_df %>% dplyr::filter(!IID %in% IDs2remove & !is.na(IID)) %>% 
  select(-data4prs$phenos2use$pheno)


#Also make a scaled df:
prsmicro_df_scale <- scale(prsmicro_df %>% select(-fgfp_id, -IID, -snpsex)) %>% as.data.frame() %>% 
  #We didnt want to scale sex and wanted to keep this as an integer, so add this variable back in
  cbind(as.factor(prsmicro_df$snpsex)) %>% rename(snpsex = "as.factor(prsmicro_df$snpsex)")

#Extract the RNT traits
RNTs <- data4prs$gwasedmts %>% dplyr::select(ends_with("RNTRes")) %>% colnames()

#Empty tibbles for output of results
#We can run this both as a multivariate and univariate regression to see how results compare
regout_RNTs_univar <- tibble()
regout_RNTs_multivar <- tibble()

for(mt in RNTs){
  #First run multivar
  m <- lm(reformulate(c(prs2use,age_sex), response = mt), data = prsmicro_df_scale)
  out <- tidy(m, conf.int = TRUE) %>% cbind(nobs(m), mt) %>% dplyr::filter(term %in% prs2use)
  regout_RNTs_multivar <- rbind(regout_RNTs_multivar, out) 
  
  #Then run it as univar
  for (prs in prs2use){
 
#Run the regression, scale it, extract coefficents of interest and bind to df
    lm2 <- lm(reformulate(c(prs, age_sex), response = mt), data = prsmicro_df_scale)
    out2 <- tidy(lm2) %>% 
      dplyr::filter(term == prs) %>% 
      cbind(nobs(lm2), mt)
  
  #Save output
  regout_RNTs_univar <- rbind(regout_RNTs_univar, out2) 
  }
}
  
#We can now save our output
regout_prsRNTs <- list(regout_RNTs_univar = regout_RNTs_univar,
                       regout_RNTs_multivar = regout_RNTs_multivar)



####################################################################
#Next we can run the observational regression of covar -> RNT MT
####################################################################
#Can just look at those with signal as highlighted above, however I think best to run for all then filter on the ones we are interested in after


#We now want to edit dataframe to remove duplicates that are present from repeat genetic measures

micropheno_df <- prsmicropheno_df %>% select(-IID, -prs2use) %>% distinct() %>%
  filter(!if_all(-fgfp_id, is.na))

#Also make a scaled df:
micropheno_df_scale <- scale(micropheno_df %>% select(-fgfp_id, -snpsex)) %>% as.data.frame() %>% 
  #We didnt want to scale sex and wanted to keep this as an integer, so add this variable back in
  cbind(as.factor(micropheno_df$snpsex)) %>% rename(snpsex = "as.factor(micropheno_df$snpsex)")


obs_est <- tibble()
for(mt in RNTs){

for(pheno in data4prs$phenos2use$pheno)  {
    lm <- lm(reformulate(c(pheno, age_sex), response = mt), data = micropheno_df_scale)
    out <- lm %>% 
      tidy() %>% 
      dplyr::filter(term == pheno) %>% 
      cbind(nobs(lm), mt)
    
    obs_est <- rbind(obs_est, out) 
}
}


#Edit this slightly so better display in plots
#Make another column so we can group the bug
obs_est <- obs_est %>%
  mutate(group = str_extract(mt, "^[^_]+"))

obs_est <- obs_est %>% 
  mutate(group = case_when(
    group == "P" ~ "Phyla",
    group == "C" ~ "Class",
    group == "O" ~ "Order",
    group == "F" ~ "Family",
    group == "G" ~ "Genus",
    group == "Div" ~ "Diversity", 
    TRUE ~ group  # Keep all other values unchanged
  ))

obs_est <- obs_est %>% 
  mutate(covar_name = case_when(
    term == "GPT_UL" ~ "Alannine_aminotransferase",
    term == "CK_UL" ~ "Creatine_kinase",
    term == "Gamma.GT_UL" ~ "Gamma_glutamyltransferase",
    term == "HDL.chol_mgdL" ~ "HDL",
    term == "Hemoglobine_gdL" ~ "Hemoglobine",
    term == "hip_circumference" ~ "Hip_circumference",
    term == "MCHC_g_dL" ~ "Mean_corp_hem_conc",
    term == "RBC_milj_mm3" ~ "Red_Blood_Cells",
    term == "Creatinine_mgdL" ~ "Creatinine",
    term == "Ureum_mgdL" ~ "Urea",
    term == "Urinezuur_mgdL" ~ "Uric_acid",
    term == "Height" ~ "Height",
    term == "Triglyceriden_mgdL" ~ "Triglycerides",
    TRUE ~ term  # Keep all other values unchanged
  ))



##########################################################################################
#Now we can save all objects to one .RData object so easier to load into the Markdown
##########################################################################################
save(regout_prsRNTs, obs_est,
     file = file.path(data.path, "data_out/data_Reg03.1.RData"))

