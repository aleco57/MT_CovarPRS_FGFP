#######################################################################
#During the conversion from bgen -> plink the naming has changed. 
#Lets make a mapping file so we know the new IDs for the plink format
######################################################################

library(dplyr)

#Data path for proj dir
data.path <- "../../../../data/PhenoPRS"

#Also read in the bgen sample file so we know the ordering of the genetic file
bgen_ids <- read.table(file.path(data.path, "greml/bgen/data.imputed.sample"),
                       header = F, as.is = TRUE, sep = " ", skip = 2)

plinkbgen_linker <- cbind(bgen_ids[,1:2], 1:2548, paste0("sample_", 0:2547))
colnames(plinkbgen_linker) <- c("bgenid1", "bgenid2", "plinkid1", "plinkid2")

#Read in linker data
linker_data <- read.table(file.path(data.path, "../FGFP_data/linker/fgfp_sampleid_linker_20200211.txt"), 
                          header = TRUE, as.is = TRUE, sep = "\t")


#Add the fgf linker to the plinkbgenlinker
plinkbgen_linker <- merge(plinkbgen_linker,
                          linker_data[,c("fgfp_id", "genetic_id")],
                          by.x = "bgenid1",
                          by.y = "genetic_id") %>% distinct()

write.table(plinkbgen_linker, file = file.path(data.path, "greml/mapping_files/bgenplinkidlink.txt"),
            quote = F, row.names = F)

