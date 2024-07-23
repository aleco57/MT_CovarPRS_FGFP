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

#Read in hughes supp file with h2 results
load(file = file.path(data.path, "data_out/h2bugs_02.2.RData"))
RNTh2_names <- paste0(RNT_h2_traits$TaxaName, "_RNTRes")

#Read in the regressions output so we can see which bug is associated with which PRS
load(file.path(data.path, "data_out/data_Reg03.1.RData"))

#Read in plinkbgen name linker file
plinkbgen_linker <- read.table(file = file.path(data.path, "greml/mapping_files/bgenplinkidlink.txt"),
                               header = T)

############################################################
# Now we can extract the phenotypes we want to run in GREML
############################################################

#Now we can extract the phenotypes of interest i.e. we are focusing on those which are heritable
#First change to df so can extract columns of interest
data4prs$gwasedmts <- as.data.frame(data4prs$gwasedmts)
phenos4greml <- data4prs$gwasedmts[,c("linker",paste0(RNT_h2_traits$TaxaName, "_RNTRes"))]

#Now merge with bgen linker so can format a phenotype file for GREML
phenos4greml <- merge(plinkbgen_linker, phenos4greml, , by.x = "fgfp_id", by.y = "linker", all.x = T)

########################
# PRS adjust phenotypes
########################
#Now we want to add on our adjusted MTs for the PRS covariate so we can calculate the heritability for these phenotypes

#First lets residualise the phenotypes we are interested in
traits2residualise <- filter(regout_prsRNTs$regout_RNTs_univar,
                             mt %in% RNTh2_names &
                               p.value < 0.05)


#Now residualise the traits
#Create a df for this regression
residual_df <- merge(data4prs$gwasedmts %>% dplyr::select(linker, all_of(traits2residualise$mt)), 
                     data4prs$pheno_covariate_prs %>% dplyr::select(fgfp_id, IID, all_of(traits2residualise$term)),
                     by.x = "linker",
                     by.y = "fgfp_id",
                     all.y = T) 

#We are only interested in those genetic_ids which are in our bgen file as have passed genetic QC
residual_df <- dplyr::filter(residual_df, 
                             IID %in% phenos4greml$bgenid1,
                             linker %in% phenos4greml$fgfp_id) %>% distinct()

#Now extract the residuals from the regression model
#Note for this we have only residualised for one PRS against one MT. 
#If multiple PRSs associated then could start looking at multiple PRSs regressed against trait

for (i in 1:nrow(traits2residualise)){
  
  bug <- traits2residualise$mt[i]
  prs <- traits2residualise$term[i]
  
  residual_df[, paste0(bug,"_", prs ,"_res")] <- lm(reformulate(prs, response = bug),
                                                   data = residual_df, na.action = na.exclude) %>% residuals()
  
}




phenos4greml <- merge(phenos4greml, 
                      residual_df %>% dplyr::select(linker, IID, ends_with(".Pt_5e.08_res")),
                      by.x = c("fgfp_id", "bgenid1"),
                      by.y = c("linker", "IID"),
                      all = T)

phenfile <- phenos4greml[,-c(1:3)]

#Finally order the ids so in the same order as plink id file
phenfile <- phenfile[order(phenfile$plinkid1),]

#Save as .phen file for GREML
write.table(phenfile,
            file.path(data.path, "greml/phenos4greml/h2mts.phen"),
            col.names = F, 
            sep = "\t",
            row.names = F,
            quote = F)

writeLines(colnames(phenfile)[-c(1:2)], file.path(data.path, "greml/phenos4greml/phennames.txt"))

#Also save as an R object as may be useful later
save(phenos4greml, file = file.path(data.path, "greml/phenos4greml/data_phenos04.2.1.RData"))


