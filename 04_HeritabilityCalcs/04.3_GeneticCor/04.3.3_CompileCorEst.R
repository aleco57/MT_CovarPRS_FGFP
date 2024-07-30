####################################################################################
#### This script will read in the GREML cor estimates and format them to a clean df
####################################################################################

library(dplyr)
library(readxl)
library(stringr)

#Data path for proj dir
data.path <- "../../../../data/PhenoPRS"



#Read in the h2 cor estimates from GREML
h2files <- list.files(file.path(data.path, "greml/greml_out/h2out_cor_dircons")) %>% str_subset(".hsq")


h2_estimates <- c()
for(h2filename in h2files){
  h2file <- read.table(file.path(data.path, "greml/greml_out/h2out_cor_dircons", h2filename), header = T, row.names = 1, sep = "\t", fill =T)
  h2file_vector <- h2file[c(9:11, 13),] %>% t() %>% c(h2filename) %>% na.omit()
  h2_estimates <- rbind(h2_estimates, h2file_vector)
}

#Make output a dataframe
h2_estimates <- as.data.frame(h2_estimates)

#Extract phenos from filenames
phenos <- sub("^h2out_|\\.hsq$", "", h2_estimates[,ncol(h2_estimates)]) %>% str_split("__", simplify = T)


h2_estimates <- cbind(h2_estimates, phenos)

#Make column name for them
colnames(h2_estimates) <- c("trait1_h2", "trait1_h2_se", "trait2_h2", "trait2_h2_se", "rg", "rg_se", "n", "filename", "pheno1", "pheno2")

#Edit pheno2 slightly to remove .hsq
h2_estimates$pheno2 <- gsub(".hsq", "", h2_estimates$pheno2)

#Save the formatted file
write.table(h2_estimates, file = file.path(data.path, "greml/greml_out/h2out_cor_dircons/clean_out.txt"), row.names = F, quote = F, sep = "\t")


#Now also make the filtered so we know which traits we want to re-run h2 esitmates for
#First generate our stable genetic correlation results
genetic_cor_filt <- filter(h2_estimates, rg_se < 0.5 & rg_se > 0)

#Load in obs_prs est
load(file.path(data.path, "data_out/data_ConsistentAssoc03.2.RData"))

#Now merge with obs and PRS results so we can further filter on consistency
obs_prs_cor <- merge(obsest_prsest,
                     genetic_cor_filt,
                     by.x = c("term", "mt"),
                     by.y = c("pheno1", "pheno2"))

#Now filter on those with consistent directionality
obs_prs_cor_filt <- obs_prs_cor[obs_prs_cor$estimate_obs*obs_prs_cor$rg > 0,]

#Save the formatted file
write.table(obs_prs_cor_filt, file = file.path(data.path, "greml/greml_out/h2out_cor_dircons/clean_out_withobsprs.txt"), row.names = F, quote = F, sep = "\t")


