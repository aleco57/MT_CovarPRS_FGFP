##################################################################################################################################
#### This script will filter to those traits we are interested in for h2 estimation and reprocess them adjusting for the env_PRSs
##################################################################################################################################

library(dplyr)
library(readxl)

#Data path for proj dir
data.path <- "../../../../data/PhenoPRS"

#Load in our data
load(file.path(data.path, "data_out/fgfpdata4prs.RData"))

#Also load in path to scratch, proj_dir and fgfp data
source("../../../parameters/scratch_path.txt")
source("../../../parameters/base_dir.R")

#Read in hughes supp file with h2 results
hughes_h2 <- read_xlsx(path = file.path(data.path, "Hughes_mtGWAS.xlsx"),
                       sheet = "Table S3",
                       skip = 4)

#Change pval cols to numeric so we can filter on these values
pvalcols <- colnames(hughes_h2)[colnames(hughes_h2) %>% endsWith("_pval")]
hughes_h2[,pvalcols] <- sapply(hughes_h2[,pvalcols], as.numeric)

#Now find out which RNT traits showed evidence for heritability
RNT_h2_traits <- filter(hughes_h2, RNT_pval < 0.05) %>% dplyr::select(TaxaName, RNT_h2, RNT_se, RNT_pval, RNT_n)
RNT_h2_traits

#Also read in the bgen sample file so we know the ordering of the genetic file
bgen_ids <- read.table(file.path(data.path, "greml/bgen/data.imputed.sample"),
                       header = F, as.is = TRUE, sep = " ", skip = 2)

#Now we can extract the phenotypes of interest
#First change to df so can extract columns of interest
data4prs$gwasedmts <- as.data.frame(data4prs$gwasedmts)
phenos <- data4prs$gwasedmts[,c("linker",paste0(RNT_h2_traits$TaxaName, "_RNTRes"))]
phenos <- merge(phenos, data4prs$pheno_covariate_prs[,c("fgfp_id", "IID")], by.x = "linker", by.y = "fgfp_id", all = T)

#Filter to only ids in the bgen file
#Use distinct as there are some repeat observations in the df
phenos <- filter(phenos, IID %in% bgen_ids$V1) %>% distinct()

#Now order it the same as the bgen file
phenos <- phenos[match(bgen_ids$V1, phenos$IID), ]
#Check the ordering worked
identical(phenos$IID, bgen_ids$V1)

#Output file for GRM  
RNT_phenofile_fgfpids <- cbind(phenos$IID, phenos$IID, phenos[,2:8])

#Also make an output file with new names after plink file conversion, these are numbers 1:2458 then sample_0:2547
RNT_phenofile_plinkid <- cbind(1:2548, paste0("sample_", 0:2547), phenos[,2:8])

#Save these pheno files for h2 of grm calculation
write.table(RNT_phenofile_fgfpids,
            file.path(data.path, "greml/phenos4greml/sigRNT_fgfpids.phen"),
            col.names = F, 
            sep = "\t",
            row.names = F,
            quote = F)

write.table(RNT_phenofile_plinkid,
            file.path(data.path, "greml/phenos4greml/sigRNT_plinkids.phen"),
            col.names = F, 
            sep = "\t",
            row.names = F,
            quote = F)

########################
# PRS adjust phenotypes
########################
#Merge the traits of interest with PRSs
phenos_prs <- merge(phenos, 
                    data4prs$pheno_covariate_prs %>% dplyr::select(ends_with(".08"), fgfp_id, IID) %>% distinct(),
                    by.x = c("linker", "IID"),
                    by.y = c("fgfp_id", "IID"),
                    all.x = T)

#See which PRS's are associated with the traits
#First generate a list of prs we are interested in
prs2use <- colnames(phenos_prs)[colnames(phenos_prs) %>% endsWith(".08") & !colnames(phenos_prs) %in% c("creatinine2.Pt_5e.08", "stool_freq.Pt_5e.08")]

#Scale PRS so all on the same scale
phenos_prs[,prs2use] <- scale(phenos_prs[,prs2use])

#Make df to store which PRSs are associated with which MT
sig_prs <- data.frame(matrix(ncol=2, nrow = 0))
for(bug in paste0(RNT_h2_traits$TaxaName, "_RNTRes")){
  #Run the lm model with all the PRS and extract residuals
  phenos_prs[, paste0(bug, "_allprs_res")] <- lm(reformulate(prs2use, response = bug),
                                                data = phenos_prs, na.action = na.exclude) %>% residuals()
  
  #Now we can filter to those PRS with evidence of an association (p < 0.05) and take the residuals of this regression
  lm_model <- lm(reformulate(prs2use, response = bug),
     data = phenos_prs, na.action = na.exclude) %>% summary()
 
  #Extract those variables which suggest evidence for 
  sigprs <- rownames(lm_model[["coefficients"]])[lm_model[["coefficients"]][,"Pr(>|t|)"] < 0.05]
  
  if(length(sigprs) > 0){
    phenos_prs[, paste0(bug, "_sigprs_res")] <- lm(reformulate(sigprs, response = bug),
                                                  data = phenos_prs, na.action = na.exclude) %>% residuals()
  }else{sigprs <- NA}
  
  sig_prs <- rbind(sig_prs, c(bug, paste(sigprs, collapse = ",")))
  
}

#Change colnames of our sig_prs df
colnames(sig_prs) <- c("bug", "sig_prs")

#Save this data for script 4.2.3
save(phenos_prs, sig_prs, file = file.path(data.path, "greml/data4script4.2.3.RData"))

#Now we can output our phenotypes ready for GREML
RNT_PRSRes_phenofile_plinkid <- cbind(1:2548, paste0("sample_", 0:2547), phenos_prs[,colnames(phenos_prs)[endsWith(colnames(phenos_prs), "_res")]])

write.table(RNT_PRSRes_phenofile_plinkid,
            file.path(data.path, "greml/phenos4greml/sigRNT_PRSRes_plinkids.phen"),
            col.names = F, 
            sep = "\t",
            row.names = F,
            quote = F)
