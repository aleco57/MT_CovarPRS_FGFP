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

write.table(plinkbgen_linker, file = file.path(data.path, "greml/mapping_files/bgenplinkidlink.txt"),
            quote = F, row.names = F)

