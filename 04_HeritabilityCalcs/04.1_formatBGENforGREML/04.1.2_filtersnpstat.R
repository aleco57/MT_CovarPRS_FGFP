#########################################
##
##  A tiny R script to identify 
##  extreme HWE allele, filter on INFO
## 	and filter on MAF from snpsum
##  qctools ouput
##
##########################################
library(dplyr)
args = commandArgs(trailingOnly=TRUE)
input = args[1]
output = args[2]

snpsum <- read.table(input, header = TRUE, as.is = TRUE, sep = " ") %>% 
filter(minor_allele_frequency > 0.05 & info > 0.9 & HW_exact_p_value > 1e-5)

writeLines(as.character(snpsum[['rsid']]), con = file(output), sep = " ")
