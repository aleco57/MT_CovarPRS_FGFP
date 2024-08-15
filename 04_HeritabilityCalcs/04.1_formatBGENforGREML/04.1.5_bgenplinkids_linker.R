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

plinkbgen_linker <- plinkbgen_linker[order(plinkbgen_linker$plinkid1),]

write.table(plinkbgen_linker, file = file.path(data.path, "greml/mapping_files/bgenplinkidlink.txt"),
            quote = F, row.names = F)

########################################################################################
#### Extract the relatedness ids we want to remove from the analysis - Dave's sample size 
########################################################################################

#Data path for proj dir
data.path <- "../../../../data"

excl <- read.table(file.path(data.path, "FGFP_data/FGFP_MT/sample_qc/Sample_Exclusion_Criteria.txt"),
                   header = TRUE, as.is = TRUE, sep = "\t")

bgen_ids <- read.table(file.path(data.path, "FGFP_data/FGFP_MT/sample_qc/data.imputed.sample"),
                       header = F, as.is = TRUE, sep = " ", skip = 2)


#Extract the relatedness ids we want to remove from the analysis
id2remove <- excl[excl$Related == 1,] %>% rownames()

plinkids2remove <- filter(plinkbgen_linker, bgenid1 %in% id2remove) %>% select(plinkid1, plinkid2)

write.table(plinkids2remove, file = file.path(data.path, "PhenoPRS/greml/plink/pruned/dave_sample/relatedness_ids2remove.txt"),
              quote = F, row.names = F, col.names = F)


plinkbgen_linker_davesample <- filter(plinkbgen_linker, !(bgenid1 %in% id2remove))

write.table(plinkbgen_linker_davesample, file = file.path(data.path, "PhenoPRS/greml/mapping_files/bgenplinkidlink_davesample.txt"),
            quote = F, row.names = F)
