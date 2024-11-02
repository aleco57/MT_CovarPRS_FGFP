#################################################################################################
#### This script will generate the phenotypes we are interested in for calculating heritability
################################################################################################

###########
# Set up
###########
library(dplyr)
library(readxl)

#Data path for proj dir
data.path <- "../../../../data/PhenoPRS"

#Also load in path to scratch, proj_dir and fgfp data
source("../../../parameters/scratch_path.txt")
source("../../../parameters/base_dir.R")

#Load in our data
load(file.path(data.path, "data_out/fgfpdata4prs.RData"))

#Read in the correlation output so we can see where we have consistent estimates
obs_prs_cor_filt <- read.table(file = file.path(data.path, "greml/greml_out/h2out_cor_dircons/clean_out_withobsprs.txt"), 
                               header = T, sep = "\t")

#Load in the data from regressions
load(file.path(data.path, "data_out/data_Reg03.1.RData"))

load(file.path(data.path, "data_out/data_ConsistentAssoc03.2.RData"))

#Read in plinkbgen name linker file
plinkbgen_linker <- read.table(file = file.path(data.path, "greml/mapping_files/bgenplinkidlink.txt"),
                               header = T)

############################################################
# Now we can extract the phenotypes we want to run in GREML
############################################################

# Edit our "data4prs[["matchedvars"]]" slightly to remove PRSs we are no longer interested in
data4reg <- data4prs[["matchedvars"]]

data4prs$gwasedmts <- as.data.frame(data4prs$gwasedmts)
phenos4greml_dircons <- data4prs$gwasedmts[,c("linker", unique(candidate_bugs$mt))]

#Now merge with bgen linker so can format a phenotype file for GREML
phenos4greml_dircons <- merge(plinkbgen_linker, phenos4greml_dircons, , by.x = "fgfp_id", by.y = "linker", all.x = T)


#Now residualise the traits
#Create a df for this regression
residual_df_dircons <- merge(data4prs$gwasedmts %>% dplyr::select(linker, all_of(candidate_bugs$mt)), 
                             data4prs$pheno_covariate_prs %>% dplyr::select(fgfp_id, IID, all_of(candidate_bugs$term_prs)),
                             by.x = "linker",
                             by.y = "fgfp_id",
                             all.y = T) 

#We are only interested in those genetic_ids which are in our bgen file as have passed genetic QC
residual_df_dircons <- dplyr::filter(residual_df_dircons, 
                                     IID %in% phenos4greml_dircons$bgenid1,
                                     linker %in% phenos4greml_dircons$fgfp_id) %>% distinct()

################## 
# Residualisation
##################


#Now extract the residuals from the regression model
#First we can do this for 1:1

for (i in 1:nrow(candidate_bugs)){
  
  bug <- candidate_bugs$mt[i]
  prs <- candidate_bugs$term_prs[i]
  
  residual_df_dircons[, paste0(bug, "_res_", prs)] <- lm(reformulate(prs, response = bug),
                                                         data = residual_df_dircons, na.action = na.exclude) %>% residuals()
  
}

#Then we can also extract the residuals for the multiple pheno residualisation of the microbial trait

multiplesignal_mt <- candidate_bugs$mt[duplicated(candidate_bugs$mt)]

for (i in 1:length(multiplesignal_mt)){
  
  bug <- multiplesignal_mt[i]
  prs <- candidate_bugs[candidate_bugs$mt == bug, "term_prs"] 
  prs_name <- paste(prs, collapse = "__")
  
  residual_df_dircons[, paste0(bug, "_res_", prs_name)] <- lm(reformulate(prs, response = bug),
                                                              data = residual_df_dircons, na.action = na.exclude) %>% residuals()
  
}


phenos4greml_dircons <- merge(phenos4greml_dircons[,1:5], 
                              residual_df_dircons[,!colnames(residual_df_dircons) %in% unique(candidate_bugs$term_prs)],
                              by.x = c("fgfp_id", "bgenid1"),
                              by.y = c("linker", "IID"),
                              all = T)

phenfile2 <- phenos4greml_dircons[,-c(1:3)]

#Finally order the ids so in the same order as plink id file
phenfile2 <- phenfile2[order(phenfile2$plinkid1),]

#Save as .phen file for GREML
write.table(phenfile2,
            file.path(data.path, "greml/phenos4greml/dircons/dircons.phen"),
            col.names = F, 
            sep = "\t",
            row.names = F,
            quote = F)

writeLines(colnames(phenfile2)[-c(1:2)], file.path(data.path, "greml/phenos4greml/dircons/dircons_phenos.txt"))

#Also save as an R object as may be useful later
save(phenos4greml_dircons, file = file.path(data.path, "greml/phenos4greml/dircons/data_phenosdircons_04.2.1.RData"))




#We are no longer using this code below











################################################################### 
# Now we can repeat the analysis for the heritable bug associations
###################################################################

#Load in which bugs are heritable
load(file.path(data.path, "data_out/h2bugs_02.2.RData"))

phenos4greml <- data4prs$gwasedmts[,c("linker", paste0(RNT_h2_traits$TaxaName, "_RNTRes"))]

#Now merge with bgen linker so can format a phenotype file for GREML
phenos4greml <- merge(plinkbgen_linker, phenos4greml, by.x = "fgfp_id", by.y = "linker", all.x = T)



########################
# PRS adjust phenotypes
########################
#Merge the traits of interest with PRSs
phenos_prs <- merge(phenos4greml, 
                    data4prs$pheno_covariate_prs %>% select(ends_with(".08"), fgfp_id, IID) %>% distinct(),
                    by.x = c("fgfp_id", "bgenid1"),
                    by.y = c("fgfp_id", "IID"),
                    all.x = T) %>% distinct()

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

#Remove columns in phenos_prs that are end with 5e.08
phenfile <- phenos_prs %>% select(-ends_with(".08"))

phenfile <- phenfile[,-c(1:3)]

#Order by plinkid1
phenfile <- phenfile[order(phenfile$plinkid1),]

#Save as .phen file for GREML
write.table(phenfile,
            file.path(data.path, "greml/phenos4greml/h2bugs/h2bugs.phen"),
            col.names = F, 
            sep = "\t",
            row.names = F,
            quote = F)

writeLines(colnames(phenfile)[-c(1:2)], file.path(data.path, "greml/phenos4greml/h2bugs/h2bugs.txt"))



