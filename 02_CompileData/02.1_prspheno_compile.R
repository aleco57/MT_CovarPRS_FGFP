#Compiles data to one file (PRS/Phenos) so can be easily used in subsequent analyses
#####################################################################################


#This script can run on RStudio
#setwd to where the script is located

#Library
library(dplyr)
library(data.table)

#Data path for proj dir
data.path <- "../../../data/PhenoPRS"

#Also load in path to scratch, proj_dir and fgfp data
source("../../parameters/scratch_path.txt")
source("../../parameters/base_dir.R")
source(file.path(pdir, "scripts/parameters/parafile.txt"))

#Load in Dave's processed fgfp data
#Load in dave processed data
load(file.path(pdir2, "data/mydata_update.Rdata"))


#Also load in the GWAS trait table we were looking at
gwas_table <- read.csv(file = file.path(data.path, "ProcessedGWASTable_OG_GC.csv"), header = T)


phenos4prs <- mydata[["pheno_data"]][["pheno_data"]][c("age", 
                                                       "BMI",
                                                       "stoole_score",
                                                       "Hemoglobine_gdL", 
                                                       "RBC_milj_mm3", 
                                                       "Creatinine_mgdL",
                                                       "e.GFR",
                                                       "HDL.chol_mgdL",
                                                       "Urinezuur_mgdL",
                                                       "Triglyceriden_mgdL",
                                                       "fruits_average_consumption_last_week",
                                                       "fruits_consumption_frequency_last_week",
                                                       "fruits_last_consumption_before_sampling",
                                                       "meat_average_consumption_last_week",
                                                       "meat_consumption_frequency_last_week",
                                                       "meat_last_consumption_before_sampling",
                                                       "coffee_average_consumption_last_week",
                                                       "coffee_consumption_frequency_last_week",
                                                       "coffee_last_consumption_before_sampling",
                                                       "alcohol_average_consumption_last_week",
                                                       "alcohol_consumption_frequency_last_week",
                                                       "alcohol_last_consumption_before_sampling",
                                                       "sugary_soda_average_consumption_last_week",
                                                       "sugary_soda_consumption_frequency_last_week",
                                                       "sugary_soda_last_consumption_before_sampling",
                                                       "sleeping_hours_per_day")]


#Read in covariates for regression
PCs <- read.table(file.path(bdir, genotypePCs), head=T, sep = "\t")

#Edit the name slightly - seemed to be missing a 0 in its linker name so couldnt merge
#sample_name_edit2 is presumed correct as no missingness when merge with this name
batch_data <- fread(file.path(bdir, microbiome_batch_data), header = TRUE, sep = "\t")
batch_data$sample_name_edit <- gsub( "\\.", "_", batch_data$SAMPLE_NAME)
batch_data$sample_name_edit2 <- gsub( "\\.", ".0", batch_data$SAMPLE_NAME)


### Load in linker File
#Note: clincal_data_available == 1 does not reduce the n down any further
linker_data <- read.table(file.path(bdir, linker), 
                          header = TRUE, as.is = TRUE, sep = "\t")
complete_data <- dplyr::filter(linker_data, micobiome_data_available == 1 & genotypedata_good_for_gwas == 1
                               & clinical_data_available == 1)

#Make covariate matrix with the outcomes we are interested in
cov_matrix <- merge(complete_data[,1:4], phenos4prs, by.x="fgfp_id", by.y = 0, all = T) %>%
  merge(PCs[,c(2,4:14)], by.x="fgfp_id", by.y = "ID_2", all = T) %>%
  merge(batch_data[,c(2:7, 11)], by.x="fgfp_id", by.y = "sample_name_edit2", all=T)

#Remove duplicated data
cov_matrix <- cov_matrix[!duplicated(cov_matrix),]

#Change date variables into years
cov_matrix$Aliquote.date_year <- strsplit(cov_matrix$Aliquote.date, "/") %>% sapply(FUN = function(x) x[[3]]) %>% as.numeric()
cov_matrix$Date.of.extraction_year <- strsplit(cov_matrix$Date.of.extraction, "/") %>% sapply(FUN = function(x) x[[3]]) %>% as.numeric()


#85% were cut and so fill in missing with cut
cov_matrix <- cov_matrix %>%
  mutate(drilled.cut = ifelse(drilled.cut == "", "cut", drilled.cut))

gwastraits <- gwas_table$file_name

#Read in our prs scores for each trait
prs <- list()
for(trait in gwastraits){
  prs[[trait]] <- read.table(file = file.path(data.path, "prs_out", trait, "prsice_out", paste0(trait, ".all_score")), header = T)
             }

#Now generate a table from this list
prs_df <- do.call(cbind, prs) %>% dplyr::select(-ends_with("FID"), -ends_with("IID"), stool_freq.IID) %>%
  rename(IID = stool_freq.IID)

pheno_covariate_prs <- merge(prs_df, linker_data, by.x = "IID", by.y = "genetic_id", all = T) %>%
  merge(cov_matrix, by = "fgfp_id", all = T) 


### Lets also load in the fgfp gwased traits
#Load in microbial data
gwased_mt_data <- fread(file.path(bdir, microbiome_GWAS_source), 
                        header = F, sep = "\t")
header <- c("linker", fread(cmd = paste("head -n 1", file.path(bdir, microbiome_GWAS_source))) %>% colnames())
colnames(gwased_mt_data) <- header



#Now lets save the clean dataframe / list into the data directory
#Some repeats of data so have used distinct() command to make sure these are removed
data4prs <- list(pheno_covariate_prs = distinct(pheno_covariate_prs), 
                 gwasedmts = gwased_mt_data,
                 rawmts = mydata[["micro_data"]])


data4prs[["matchedvars"]] <- data.frame(prs = c(paste0(gwas_table$file_name, ".Pt_5e.08"), "creatinine2.Pt_5e.08"),
                                  pheno = c("stoole_score",
                                            "Hemoglobine_gdL",
                                            "RBC_milj_mm3",
                                            "Urinezuur_mgdL",
                                            "HDL.chol_mgdL",
                                            "BMI",
                                            "Creatinine_mgdL",
                                            "Triglyceriden_mgdL",
                                            "sugary_soda_average_consumption_last_week",
                                            "meat_average_consumption_last_week",
                                            "coffee_average_consumption_last_week",
                                            "sleeping_hours_per_day",
                                            "alcohol_average_consumption_last_week",
                                            "fruits_average_consumption_last_week",
                                            "Creatinine_mgdL",
                                            "e.GFR"))

save(data4prs, file = file.path(data.path, "data_out/fgfpdata4prs.RData"))
