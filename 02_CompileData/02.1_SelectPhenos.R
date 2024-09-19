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
load(file.path(pdir2, "data/mydata.Rdata"))


#Also load in the GWAS trait table we were looking at
#gwas_table <- read.csv(file = file.path(data.path, "ProcessedGWASTable_OG_GC.csv"), header = T)


#Make a search term to find variables of interest
grep("cir", colnames(mydata[["pheno_data"]]), value = TRUE, ignore.case = TRUE)

phenfromdave <- mydata[["pheno_data"]][c("stoole_score",
                                         "RBC_milj_mm3",
                                         "Urinezuur_mgdL", #Uric acid
                                         "Hemoglobine_gdL",
                                         "e.GFR",
                                         "BMI",
                                         "age", 
                                         "HDL.chol_mgdL",
                                         "J01CA04",
                                         "Height",
                                         "Triglyceriden_mgdL",
                                         "fruits_last_consumption_before_sampling",
                                         "coffee_last_consumption_before_sampling",
                                         "beer_comsuption",
                                         "alcohol_average_consumption_last_week",
                                         "meat_average_consumption_last_week",
                                         "Gluc_nuchter_mg.dL", "Insuline_nuchter_mUL", #Work out HOMA-IR
                                         "sleeping_hours_per_day",
                                         "soy_products_average_consumption_last_week",
                                         "GPT_UL", #This is ALT
                                         "CK_UL", #Creatine Kinase
                                         "G03CA04_CC06",
                                         "Creatinine_mgdL",
                                         "Gamma.GT_UL",
                                         "MCHC_g_dL",
                                         "A07EC02",
                                         "Ureum_mgdL", #Urea
                                         "FOODSUPPLEMENT_Q10",
                                         "L04AX01",
                                         "N06AX16",
                                         "R06AX28",
                                         "PARTICIPANTS_ALLERGY_food", "PARTICIPANTS_ALLERGY_food.1",
                                         "FOODSUPPLEMENT_Mg")]

#The following variables are missing so we will search them from the excel sheet provided:

#Cigarettes average per day missing
#Hip circumference missing
#Chronic Fatigue Syndrome missing
#Irritable Bowel Syndrome missing
#Soking status missing
#Ulcerative Colitis missing

#First load the excel spreadsheet
excel_pheno <- readxl::read_xlsx(file.path(bdir, "data/fgf/phenotypic/clinic_data/released/2024-07-29/data/fgfp_20240507/metadata.fgfp.xlsx")) %>%
  select(vdp_ids, average_cigarettes_consumption_per_day_last3months, hip_circumference, 
         chronic_fatigue_syndrome, ibs, have_smoked_last3month, had_ulcerative_colitis)

#To merge these we also need the linker file that dave has made
### Load in linker File
#Note: clincal_data_available == 1 does not reduce the n down any further
linker_data <- read.table(file.path(bdir, linker), 
                          header = TRUE, as.is = TRUE, sep = "\t")
complete_data <- dplyr::filter(linker_data, micobiome_data_available == 1 & genotypedata_good_for_gwas == 1
                               & clinical_data_available == 1)

allphens <- left_join(phenfromdave %>% mutate(fgfp_id = rownames(phenfromdave)), 
                      unique(linker_data[,c("fgfp_id","labid_16S")]), by = "fgfp_id") %>%
  left_join(excel_pheno %>% mutate(edit_vdp_ids = gsub("\\.(?=[0-9])", ".0", vdp_ids, perl = TRUE)), 
            by = c("fgfp_id" = "edit_vdp_ids"))

#We also want to work out HOMA-IR as this is not present in the df
allphens <- allphens %>%
  mutate(HOMA_IR = (Gluc_nuchter_mg.dL * Insuline_nuchter_mUL) / 405) %>% select(-Gluc_nuchter_mg.dL, -Insuline_nuchter_mUL) %>%
  select(fgfp_id, vdp_ids, labid_16S, everything())

#Now we want to save the phenotype data for later use
write.csv(distinct(allphens), file = file.path(data.path, "data_out/fgfp_phenotypes4prs.csv"), row.names = F, quote = F)
