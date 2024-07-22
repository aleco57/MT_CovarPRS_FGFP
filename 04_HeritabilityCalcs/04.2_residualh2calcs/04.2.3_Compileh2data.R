####################################################################################
#### This script will read in the GREML h2 estimates and format them to a clean df
####################################################################################

library(dplyr)
library(readxl)
library(stringr)

#Data path for proj dir
data.path <- "../../../../data/PhenoPRS"

#Load data from previous script
load(file.path(data.path, "greml/phenos4greml/data_phenos04.2.1.RData"))

#Read in the h2 estimates from GREML
#First read in the unadjusted h2
h2files <- list.files(file.path(data.path, "greml/greml_out/h2out")) %>% str_subset(".hsq")

h2_estimates <- c()
for(h2filename in h2files){
  h2file <- read.table(file.path(data.path, "greml/greml_out/h2out", h2filename), header = T, row.names = 1, sep = "\t", fill =T)
  h2_estimates <- rbind(h2_estimates, c(h2filename, h2file["V(G)/Vp","Variance"], h2file["V(G)/Vp","SE"], h2file["Pval","Variance"]))
}
#Make output a dataframe
h2_estimates <- as.data.frame(h2_estimates)

#Match these to the correct phenotypes
colnames(h2_estimates) <- c("file", "h2", "se", "pval")
h2_estimates$trait <- gsub("h2out_", "", h2_estimates$file)
h2_estimates$trait <- gsub(".hsq", "", h2_estimates$trait)

#Save the formatted file
write.table(h2_estimates, file = file.path(data.path, "greml/greml_out/h2out/cleanh2est.txt"), row.names = F, quote = F, sep = "\t")
