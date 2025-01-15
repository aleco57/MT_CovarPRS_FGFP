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

#Make a list of unique prs traits
prstraits <- unique(candidate_bugs$prs)


#First make df to store results
bug_sumstats <- list()

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
    merged_gwasbetas[[covar]][[bug]] <- merge(bug_sumstats[[covar]][[bug]],
                      snps_df[[covar]],
                      by.x = "snpid", by.y = "fgfp_linker")
  }
}

save(bug_sumstats, merged_gwasbetas, file = file.path(data.path, "data_out/filtered_bug_sumstats_05.1.RData"))


##########################################################################################################
##########################################################################################################
##########################################################################################################

### Now we want to repeat but looking at the mibiogen snps
mbgbugs <- c("simpson", #Div_NumGen
             "invsimpson", #Div_Chao
             "family.Lachnospiraceae.id.1987",
             "genus.Barnesiella.id.944",
             "genus.Roseburia.id.2012", 
             "family.Ruminococcaceae.id.2050", #G_Sporobacter
             "family.Ruminococcaceae.id.2050", #G_unclassified
             NA, #G_Unclassified__K_Bacteria
             "order.Clostridiales.id.1863",
             "phylum.Firmicutes.id.1672", #G_Unclassified
             "family.Erysipelotrichaceae.id.2149",
             "genus.Holdemanella.id.11393",
             "genus.Ruminococcus1.id.11373",
             "family.Erysipelotrichaceae.id.2149", #G_unclassified
             "family.Ruminococcaceae.id.2050",
             "genus.Oscillibacter.id.2063",
             "genus.Bacteroides.id.918",
             "family.Rhodospirillaceae.id.2717", #G_Aestuariispira
             "family.Lachnospiraceae.id.1987", #G_Coprococcus
             "shannon",
             "family.Porphyromonadaceae.id.943")

matched_bugs <- data.frame(fgfpbugs = unique(candidate_bugs$mt), mbgbugs = mbgbugs)

#Remove the NA
matched_bugs <- matched_bugs[!is.na(matched_bugs$mbgbugs),]

#Path to where to nofgf sumstats are
nofgf.path <- "../../../data/mibiogen_data/mibiogen_betafiles/nofgf_withMAF"

#We again want to repeat the code above to find matching
bug_sumstats_mbg <- list()

for(prs in prstraits){
  #Extract rsids of interest for prs
  snps <- snps_df[[prs]]$rsid
  
  #Extract bugs associated with this prs
  bugs <- candidate_bugs[candidate_bugs$prs==prs, "mt"]
  
  #Remove those with no mbg sumstats
  mbg_bug <- matched_bugs$mbgbugs[matched_bugs$fgfpbugs %in% bugs]
  
  
  for(bug in mbg_bug){
    #Read in the sumstats
    
    #Slightly different path for div metrics
    if(bug %in% c("simpson", "invsimpson", "shannon")){
      mbgbugss <- file.path(nofgf.path, paste0(bug, "_EURnofgf.summary.txt.gz"))} else{
    mbgbugss <- file.path(nofgf.path, paste0(bug, ".summaryMAF.txt.gz"))
      }
    if(file.exists(mbgbugss)){
    bug_sumstats_mbg[[prs]][[bug]] <- data.table::fread(file = mbgbugss) %>% 
      filter(rsID %in% snps)
    }else{print(paste0("No sumstats for ", bug))}
    }
  }


#Also merge with relevant bug sumstats so can compare beta coefficents
merged_gwasbetas_mbg <- list()
for(covar in names(bug_sumstats_mbg)){
  for(bug in names(bug_sumstats_mbg[[covar]])){
    merged_gwasbetas_mbg[[covar]][[bug]] <- merge(bug_sumstats_mbg[[covar]][[bug]],
                                              snps_df[[covar]],
                                              by.x = "rsID", by.y = "rsid")
  }
}

save(bug_sumstats_mbg, merged_gwasbetas_mbg, file = file.path(data.path, "data_out/filtered_bug_sumstatsMBG_05.1.RData"))

