#### This script will read in the candidate bugs in and the SNPs for the covars and extract the relavent SNPs from Dave's GWAS files
# This script must be run on bp

#Library
library("dplyr")
library("data.table")

#Data path for proj dir
data.path <- "../../../data/PhenoPRS"

#Load in the instrument data for the covars
load(file.path(data.path, "phenosumstats/filtered_sumstats_list.RData"))

#Load in the candidate bugs
load(file.path(data.path, "data_out/data_ConsistentAssoc03.2.RData"))


#Rename snps file with the new naming
#Specify a function for cleaning directory names
clean_directory_name <- function(x) {
  x %>%
    stringr::str_replace_all("[^[:alnum:]]+", "_") %>% 
    stringr::str_replace_all("_+", "_")   
}


names(snps_df) <- clean_directory_name(names(snps_df))


#First make df to store results
bug_sumstats <- list()
prstraits <- unique(candidate_bugs$prs)

for(prs in prstraits){
snps <- snps_df[[prs]]$fgfp_linker 

bugs <- candidate_bugs[candidate_bugs$prs==prs, "mt"]
bugs <- gsub("Res", "", bugs)
for(bug in bugs){
bug_sumstats[[prs]][[bug]] <- data.table::fread(file = file.path(data.path, "davesRNTsumstats", paste0(bug, "_allchr.txt.gz"))) %>% 
  filter(snpid %in% snps)
}
}

#Also merge with relevant bus sumstats so can compare beta coefficents
merged_gwasbetas <- list()
for(covar in names(bug_sumstats)){
  for(bug in names(bug_sumstats[[covar]])){
    merged_gwasbetas[[covar]][[bug]] <- merge(bug_sumstats[[covar]][[bug]][,c("snpid", "frequentist_add_beta_1", "frequentist_add_se_1", "frequentist_add_pvalue")],
                      snps_df[[covar]],
                      by.x = "snpid", by.y = "fgfp_linker")
  }
}

save(bug_sumstats, merged_gwasbetas, file = file.path(data.path, "data_out/filtered_bug_sumstats_05.1.RData"))









