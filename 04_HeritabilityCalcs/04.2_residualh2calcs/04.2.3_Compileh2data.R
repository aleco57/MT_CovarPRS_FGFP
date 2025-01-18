####################################################################################
#### This script will read in the GREML h2 estimates and format them to a clean df
####################################################################################

library(dplyr)
library(readxl)
library(stringr)

#Data path for proj dir
data.path <- "../../../../data/PhenoPRS"

i <- "dircons"
  #Read in the h2 estimates from GREML
  h2files <- list.files(file.path(data.path, "greml/greml_out/h2_out", i)) %>% str_subset(".hsq")
  
  h2_estimates <- c()
  for(h2filename in h2files){
    h2file <- read.table(file.path(data.path, "greml/greml_out/h2_out",i, h2filename), header = T, row.names = 1, sep = "\t", fill =T)
    h2_estimates <- rbind(h2_estimates, c(h2filename, h2file["V(G)/Vp","Variance"], h2file["V(G)/Vp","SE"], h2file["Pval","Variance"]))
  }
  #Make output a dataframe
  h2_estimates <- as.data.frame(h2_estimates)
  colnames(h2_estimates) <- c("file", "h2", "se", "pval")
  
  #Order the file column in order of number
  h2_estimates$file <- gsub("h2out_", "", h2_estimates$file)
  h2_estimates$file <- gsub(".hsq", "", h2_estimates$file)
  h2_estimates$file <- as.numeric(h2_estimates$file)
  h2_estimates <- h2_estimates[order(h2_estimates$file),]
  
  #Match these to the correct phenotypes
  
  #Load in the phenotype names
  name_vector <- readLines(file.path(data.path, "greml/phenos4greml/dircons/dircons_phenos.txt"))
  h2_estimates$trait <- name_vector
  
  #Save the formatted file
  write.table(h2_estimates, file = file.path(data.path, "greml/greml_out/h2_out/", i, "cleanh2est.txt"), row.names = F, quote = F, sep = "\t")
  
