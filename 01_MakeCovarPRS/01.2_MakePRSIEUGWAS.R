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