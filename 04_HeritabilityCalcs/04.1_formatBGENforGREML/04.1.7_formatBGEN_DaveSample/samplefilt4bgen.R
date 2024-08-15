#Script to see which individuals have been filtered out of bgen file
source("../parameters/base_dir.R")

library("dplyr")

excl <- read.table(file.path(pdir, "data/FGFP_MT/sample_qc/Sample_Exclusion_Criteria.txt"),
           header = TRUE, as.is = TRUE, sep = "\t")

bgen_ids <- read.table(file.path(pdir, "data/FGFP_MT/sample_qc/data.imputed.sample"),
                       header = F, as.is = TRUE, sep = " ", skip = 2)

#Repeat samples
###############

rownames(excl) %>% head()
#The rownames correspond to the genetic_id which are all unique values
excl$id %>% head()
#This is the FGFP tag and there are some repeats of this

unique(bgen_ids$V1) %>% length() == length(bgen_ids$V1)
#This statement is true

#Are the two id columns the same?
identical(bgen_ids$V1, bgen_ids$V2)
#yes they are


#What are the extra identification tags in exclusion file?
identical(excl[,1], excl[,2])
#This is false meaning there are some ids which do not match
excl[excl$id != excl$id2,]
#Seems to be a two rows which are labelled "RVDP" instead of VDP

## Now lets generate a list of names for each QC step where samples were filtered out
idsfromexcl <- list()
for(x in colnames(excl)[3:ncol(excl)]){
  idsfromexcl[[x]] <-  rownames(excl)[excl[,x] == 1]  
}


#Now we can see what has been removed from the bgen files
lapply(idsfromexcl, function(x) bgen_ids$V1 %in% x %>% sum(na.rm=T))

var_notfilt <- lapply(idsfromexcl, function(x) bgen_ids$V1 %in% x %>% sum(na.rm=T)) > 0
names(idsfromexcl)[var_notfilt]
# [1] "duplicates"       "FscalePC"         "Related"          "nomicrobiomedata"

#It appears that some of the duplicates have been filtered out, there are 64 ids present in the bgen file
# and there are 80 that are flagged as duplicates
# Is it likely that the 16 removed were done so due to other filters?

dupesremoved <- idsfromexcl$duplicates[!idsfromexcl$duplicates %in% bgen_ids$V1]

excl[rownames(excl) %in% dupesremoved,-c(1:2)] %>% rowSums() > 1
#We can see that all these values are true and so the duplicates that have been removed are due to other qc factors

### We see the same pattern for the relatedness where 7/262 have been removed: is this because of other qc?
relatremoved <- idsfromexcl$Related[!idsfromexcl$Related %in% bgen_ids$V1]
excl[rownames(excl) %in% relatremoved,-c(1:2)] %>% rowSums() > 1

#Yes

#Are all of the replicated ids also contained within the related variable?
idsfromexcl$duplicates %in% idsfromexcl$Related %>% sum(na.rm = T) 
# [1] 35

#No, there are only 35 / 80 matches

#Finally FScalePC does not need to removed, this is explained in Dave's QC file.

################################
# Mismatch in numbers !!!
################################

#This mismatch comes as relatedness has not been removed
# 262 - 7 (255) is the number of related individuals still left in the analysis.
# FGFP paper states 2293 were used for the analysis
# 2548 genetic_ids in the bgen files
# 2548 - 255 = 2293 !!
