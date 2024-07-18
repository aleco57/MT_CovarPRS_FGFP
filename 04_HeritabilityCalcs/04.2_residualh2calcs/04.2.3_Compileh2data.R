####################################################################################
#### This script will read in the GREML h2 estimates and format them to a clean df
####################################################################################

library(dplyr)
library(readxl)

#Data path for proj dir
data.path <- "../../../../data/PhenoPRS"

#Load data from previous script
load(file.path(data.path, "greml/data4script4.2.3.RData"))

#Read in the h2 estimates from GREML
#First read in the unadjusted h2
h2_estimates <- c()
for(h2filename in paste0("h2out_pheno", 1:7, ".hsq")){
  h2file <- read.table(file.path(data.path, "greml/h2out", h2filename), header = T, row.names = 1, sep = "\t", fill =T)
  
  h2_estimates <- rbind(h2_estimates, c(h2filename, h2file["V(G)/Vp","Variance"], h2file["V(G)/Vp","SE"], h2file["Pval","Variance"]))
}

#Now do the same for adjusted ones
for(h2filename in paste0("h2outadj_pheno", 1:13, ".hsq")){
  h2file <- read.table(file.path(data.path, "greml/h2out/prsadj", h2filename), header = T, row.names = 1, sep = "\t", fill =T)
  
  h2_estimates <- rbind(h2_estimates, c(paste0("adj", h2filename), h2file["V(G)/Vp","Variance"], h2file["V(G)/Vp","SE"], h2file["Pval","Variance"]))
}

#Match these to the correct phenotypes
colnames(h2_estimates) <- c("file", "h2", "se", "pval")
h2_estimates <- as.data.frame(h2_estimates)
h2_estimates$pheno <- colnames(phenos_prs)[c(3:9, 25:37)]

#Calculate 95%CI for the h2 esitmates: note this is likely to be large as n is small and < suggest samples size from greml
h2_estimates <- merge(h2_estimates, sig_prs, by.x = "pheno", by.y = "bug", all=T)

#Save the formatted file
write.table(h2_estimates, file = file.path(data.path, "greml/h2out/cleanh2est.txt"), row.names = F, quote = F, sep = "\t")
