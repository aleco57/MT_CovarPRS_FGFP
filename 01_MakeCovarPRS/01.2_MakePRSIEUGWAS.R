#### Read in the GWAS table so we can extract SNPs for our instruments using IEUGWAS

#Set up
library(dplyr)
library(ieugwasr)
library(stringr)

#Data path for proj dir
data.path <- "../../../data/PhenoPRS"

#Also load in path to scratch in case need to access files here
source("../../parameters/scratch_path.txt")

#Read in the GWAS table
gwas_ref <- read.csv(file.path(data.path, "phenogwas_map.csv"), header = T, skip = 1)

#How many are not IEU?
filter(gwas_ref, IEU.CAT == "CAT")

# We will be generating 36 covariate PRSs, so now extract these
id_vector <- filter(gwas_ref, !is.na(Identifier) & IEU.CAT == "IEU") %>% select(Identifier) %>% pull()

#We will extract top hits / load in file if have already done so
if(!file.exists(file.path(data.path, "phenosumstats/allhits.csv"))){
  hits <- tophits(id=id_vector)
  write.csv(hits, file.path(data.path, "phenosumstats/allhits.csv"), row.names = F)
  }else{
    hits <- read.csv(file.path(data.path, "phenosumstats/allhits.csv"))
  }

#split this by trait so easier to work with
split_trait <- split(hits, hits$trait)

#Lets see how many SNPs are availble in the FGFP dataset
#Load in fgfp snps
load(file.path(data.path, "all_snpswrsids.RData"))

#Make filter and merge function + df then save output
filtandmerge <- function(x){
  snps_df <- merge(data.frame(rsid = x[["rsid"]],
                              eff.allele = x[["ea"]],
                              EAF = x[["eaf"]],
                              beta = x[["beta"]], 
                              SE = x[["se"]],
                              pval = x[["p"]]),
                   merged[c("V1", "rsid")], by = "rsid", all.x = T) %>% 
    rename(fgfp_linker = V1)
}

if(!file.exists(file.path(data.path, "phenosumstats/filtered_sumstats_list.RData"))){

snps_df <- lapply(split_trait, filtandmerge)
snps_df <- lapply(snps_df, function(x) x[!is.na(x[["fgfp_linker"]]),])

#Save the output as an RData file
save(snps_df, file = file.path(data.path, "phenosumstats/filtered_sumstats_list.RData"))
} else{
  load(file.path(data.path, "phenosumstats/filtered_sumstats_list.RData"))
}

#Specify a function for cleaning directory names
clean_directory_name <- function(x) {
  x %>%
    str_replace_all("[^[:alnum:]]+", "_") %>% 
    str_replace_all("_+", "_")   
}


# Now we have these snps, we can save these as files in a new directory
for (i in 1:length(snps_df)){
  name <- names(snps_df)[i] %>% clean_directory_name()
  outputdir <- file.path(data.path, "phenosumstats",  name)
  dir.create(outputdir)
  outputfile <- file.path(outputdir, paste0(name, "_processedgwas.tsv"))
  write.table(snps_df[[i]], file = outputfile, quote = F, row.names = F)

#Also write out a file with just the snps we are interested in for FGFP
outputfile_snps <- file.path(outputdir, paste0(name, "_snps.txt"))
write.table(snps_df[[i]]["fgfp_linker"], file = outputfile_snps, col.names=F, row.names=F, quote=F)
}

#Save a file with all the directory names so we can loop round this in the next script
dirnames <- names(snps_df) %>% clean_directory_name()
write.table(dirnames, file = file.path(data.path, "phenosumstats/phenosumstats_dirnames.txt"), col.names=F, row.names=F, quote=F)


#A few mismatches in red blood cell count and triglycerides, as multiple effect alleles at the same point mutation, we can fix this
for(trait in c("Red blood cell count", "triglycerides")){
allele_df <- str_sub(snps_df[[trait]]$fgfp_linker, -3) %>% str_split("_", simplify = T)
vec <- snps_df[[trait]]$eff.allele

#Find bad matches
no_match_indices <- c()
# Loop through the vector and check for matches, storing indices where no match is found
for(i in 1:length(vec)){
  if(sum(allele_df[i,] %in% vec[i]) < 1){
    no_match_indices <- c(no_match_indices, i)  # Append index to the vector
  }
}

#Now we can remove these from our df so prsice can run
snps_df[[trait]] <- snps_df[["Red blood cell count"]][-no_match_indices,]

#Re write out
name <- clean_directory_name(trait)
outputdir <- file.path(data.path, "phenosumstats",  name)
outputfile <- file.path(outputdir, paste0(name, "_processedgwas.tsv"))
write.table(snps_df[[trait]], file = outputfile, quote = F, row.names = F)
}





################# Add in yengo_bmi snps to see if this improves R2 instead of GWAS on IEUGWAS
yengo_bmi <- read.table(file.path(data.path, "yengo_BMI_656.txt"), header = T)
yengo_bmi <- yengo_bmi %>% select(SNP, Tested_Allele, Other_Allele, Freq_Tested_Allele_in_HRS, BETA, SE, P)
colnames(yengo_bmi) <- c("rsid", "eff.allele", "Oth_al", "EAF", "beta", "SE", "pval")

merg_yengo <- merge(yengo_bmi,
                               merged[c("V1", "rsid")], by = "rsid", all.x = T) %>% 
  rename(fgfp_linker = V1)

#Now save this output so can make PRS
snps_df[["yengo_bmi"]] <- merg_yengo
save(snps_df, file = file.path(data.path, "phenosumstats/filtered_sumstats_list.RData"))

name <-  clean_directory_name("yengo_bmi")
outputdir <- file.path(data.path, "phenosumstats",  name)
dir.create(outputdir)
outputfile <- file.path(outputdir, paste0(name, "_processedgwas.tsv"))
write.table(snps_df[[name]], file = outputfile, quote = F, row.names = F)

#Also write out a file with just the snps we are interested in for FGFP
outputfile_snps <- file.path(outputdir, paste0(name, "_snps.txt"))
write.table(snps_df[[name]]["fgfp_linker"], file = outputfile_snps, col.names=F, row.names=F, quote=F)

