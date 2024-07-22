#################################################################################
#A small script for extracting those MTs which showed evidence for heritability
#The script reads in results from Dave's Supplementary excel file from GWAS paper
#################################################################################
#Data path for proj dir
data.path <- "../../../data/PhenoPRS"
#Library
library(dplyr)
library(readxl)
#Read in hughes supp file with h2 results
hughes_h2 <- read_xlsx(path = file.path(data.path, "Hughes_mtGWAS.xlsx"),
                       sheet = "Table S3",
                       skip = 4)

#Change pval cols to numeric so we can filter on these values
pvalcols <- colnames(hughes_h2)[colnames(hughes_h2) %>% endsWith("_pval")]
hughes_h2[,pvalcols] <- sapply(hughes_h2[,pvalcols], as.numeric)

#Now find out which RNT traits showed evidence for heritability
RNT_h2_traits <- filter(hughes_h2, RNT_pval < 0.05) %>% dplyr::select(TaxaName, RNT_h2, RNT_se, RNT_pval, RNT_n)

save(RNT_h2_traits, file = file.path(data.path, "data_out/h2bugs_02.2.RData"))
