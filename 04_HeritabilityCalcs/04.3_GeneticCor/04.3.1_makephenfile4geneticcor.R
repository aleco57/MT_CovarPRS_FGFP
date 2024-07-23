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

#Read in hughes supp file with h2 results
load(file = file.path(data.path, "data_out/h2bugs_02.2.RData"))
RNTh2_names <- paste0(RNT_h2_traits$TaxaName, "_RNTRes")

#Extract vector of the traits so that phenotypes are in the correct order for .phen GREML file
phenovector <- lapply(1:nrow(obs_est), function(i) obs_est[i, c("term", "mt")]) %>% unlist() %>% unname()

#First lets make a df with all our variables of interest, this will be our confounders, cont PRSs and MTs
phenos4cor <- merge(data4prs$pheno_covariate_prs[,c(obs_est$term, "fgfp_id", "IID")],
                          dplyr::select(data4prs$gwasedmts, linker, obs_est$mt),
                          by.x = "fgfp_id",
                          by.y = "linker",
                          all = T) %>% 
  #Then filter to the only individuals in the bgen file
  dplyr::filter(IID %in% plinkbgen_linker$bgenid1,
                fgfp_id %in% plinkbgen_linker$fgfp_id) %>% 
  distinct()

#Now merge to make our file, must merge on both ids as some repeated measures in genetic data
phenfile <- merge(plinkbgen_linker, phenos4cor, , by.x = c("fgfp_id", "bgenid1"), by.y = c("fgfp_id", "IID")) %>%
#Sort on plinkid
  arrange(plinkid1) %>%
#Then extract columns of interest
  select(plinkid1, plinkid2, phenovector)

#Save as .phen file for GREML
write.table(phenfile,
            file.path(data.path, "greml/phenos4greml/bivariateGREML/phen4cor.phen"),
            col.names = F, 
            sep = "\t",
            row.names = F,
            quote = F)

#Also make a note of the headers for GREML
writeLines(colnames(phenfile)[-c(1:2)], file.path(data.path, "greml/phenos4greml/bivariateGREML/phennames.txt"))
