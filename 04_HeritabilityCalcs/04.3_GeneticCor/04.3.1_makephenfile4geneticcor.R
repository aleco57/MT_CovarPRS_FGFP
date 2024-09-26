##### Generting .phen file for bivarate GREML

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

#Load in regression output so we can decide which traits to run bivariate regression for
load(file.path(data.path, "data_out/data_Reg03.1.RData"))

#Read in plinkbgen name linker file
plinkbgen_linker <- read.table(file = file.path(data.path, "greml/mapping_files/bgenplinkidlink.txt"),
                               header = T)

#Read in our dir_cons associations
load(file = file.path(data.path, "data_out/data_ConsistentAssoc03.2.RData"))

#############################################################################
# We want to generate .phen file dir cons analysis
# Have now removed code which focused on h2 phenotypes from FGFP GWAS results
############################################################################

#Extract vector of the traits so that phenotypes are in the correct order for .phen GREML file
phenovector <- lapply(1:nrow(candidate_bugs), function(i) candidate_bugs[i, c("term", "mt")]) %>% unlist() %>% unname()

#First lets make a df with all our variables of interest, this will be our cont PRSs and MTs
phenos4cor <- merge(data4prs$pheno_covariate_prs[,c(unique(candidate_bugs$term), "fgfp_id", "IID")],
                    dplyr::select(data4prs$gwasedmts, linker, unique(candidate_bugs$mt)),
                    by.x = "fgfp_id",
                    by.y = "linker",
                    all = T) %>% 
  #Then filter to the only individuals in the bgen file
  dplyr::filter(IID %in% plinkbgen_linker$bgenid1,
                fgfp_id %in% plinkbgen_linker$fgfp_id) %>% 
  distinct()

#Lets scale our pheno variables to see if aids for convergence of GREML
ncols <- length(unique(candidate_bugs$term))+1
phenos4cor[,2:ncols] <- scale(phenos4cor[,2:ncols])

#Now merge to make our file, must merge on both ids as some repeated measures in genetic data
phenfile <- merge(plinkbgen_linker, phenos4cor, , by.x = c("fgfp_id", "bgenid1"), by.y = c("fgfp_id", "IID")) %>%
  #Sort on plinkid
  arrange(plinkid1)
  
#Then extract columns of interest
phenfile <- phenfile[,c("plinkid1","plinkid2",phenovector)]

#Save as .phen file for GREML
write.table(phenfile,
            file.path(data.path, "greml/phenos4greml/bivariateGREML/phen4cor_dircons.phen"),
            col.names = F, 
            sep = "\t",
            row.names = F,
            quote = F)

#Also make a note of the headers for GREML
writeLines(phenovector, file.path(data.path, "greml/phenos4greml/bivariateGREML/phennames_dircons.txt"))



