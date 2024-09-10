### Try and fix the bad PRS generated

#Set up
library(dplyr)
library(ieugwasr)

#Data path for proj dir
data.path <- "../../../data/PhenoPRS"

#Also load in path to scratch in case need to access files here
source("../../parameters/scratch_path.txt")

#Lets load in the GWAS table we have generated
load(file.path(data.path, "gwastable_gwascat.RData"))

#make vector with index of the gwas_table which didn't work
badprsvector <- c(5,10,12,15)

#Lets look in OpenGWAS instead to generate these PRS
uric_acid_ieu <- "ebi-a-GCST90018977"
alcohol_ieu <- "ukb-b-5779"
fruit_ieu <- "ebi-a-GCST90096912"
creatinine_ieu <- "ebi-a-GCST90103634"
hits <- tophits(id=c(uric_acid_ieu, alcohol_ieu, fruit_ieu))
hits2 <- tophits(id=creatinine_ieu)

#These hits are also already clumped

#Lets split these into two dataframes and compatible with the prev output
split_trait <- split(hits, hits$trait)

gwas_table <- read.csv(file = file.path(data.path, "PhenoGWASSum.csv"), header = T)
#Lets also add to the table for how many SNPs after clump
numsnps <- function(x) {
  read.table(file = file.path(data.path, "prs_out", x, paste0(x, "_snps.txt")),
             header = F) %>% nrow()
}

gwas_table$numclumpedsnps <- sapply(gwas_table$file_name, numsnps) %>% unname()

save(gwas_table, file = file.path(data.path, "gwastable_gwascat.RData"))

#Load in fgfp snps
load(file.path(data.path, "all_snpswrsids.RData"))

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


snps_df <- lapply(split_trait, filtandmerge)
snps_df <- lapply(snps_df, function(x) x[!is.na(x[["fgfp_linker"]]),])

snps_df2 <- filtandmerge(hits2)
snps_df2 <- snps_df2[!is.na(snps_df2[["fgfp_linker"]]),]


#Removing fish and plant and generating a new alcohol GWAS from IEUGWAS
#We can not edit the GWAS table and make a new one
gwas_table <- gwas_table[-c(badprsvector),]

info <- gwasinfo(c(uric_acid_ieu, alcohol_ieu, fruit_ieu))
info2 <- gwasinfo(creatinine_ieu)

gwas_table[4,"ieu_tag"] <- uric_acid_ieu
gwas_table$numclumpedsnps[4] <- nrow(snps_df$`Serum uric acid levels`)

#Make a new row for the new alcohol variable generated
gwas_table[13,] <- c("alcohol", "alcohol", "alcohol", NA, info$author[2], info$sample_size[2], NA, nrow(snps_df$`Alcohol intake frequency.`), alcohol_ieu)
gwas_table[14,] <- c("fruit", "fruit", "fruit", NA, info$author[3], info$sample_size[3], NA, nrow(snps_df$`Fresh fruit consumption`), fruit_ieu)
gwas_table[15,] <- c("creatinine3", "creatinine3", "creatinine3", NA, info2$author, info2$sample_size, NA, nrow(snps_df2), creatinine_ieu)

#Output this new df
write.csv(gwas_table, file = file.path(data.path, "ProcessedGWASTable_OG_GC.csv"), row.names = F, quote = F)



#Now lets output our SNP df

for (i in c(4,13, 14)){
if(i == 4){name = names(snps_df)[3]} else if (i == 13){name = names(snps_df)[1]}else if (i == 14){name = names(snps_df)[2]}

outputdir <- file.path(data.path, "prs_out",  gwas_table$file_name[i])
outputfile <- file.path(outputdir, paste0(gwas_table$file_name[i], "_processedgwas.tsv"))
write.table(snps_df[[name]], file = outputfile, quote = F, row.names = F)

#Also write out a file with just the snps we are interested in for FGFP
outputfile_snps <- file.path(outputdir, paste0(gwas_table$file_name[i], "_snps.txt"))
write.table(snps_df[[name]]["fgfp_linker"], file = outputfile_snps, col.names=F, row.names=F, quote=F)
}

i =15
outputdir <- file.path(data.path, "prs_out",  gwas_table$file_name[i])
dir.create(outputdir)
outputfile <- file.path(outputdir, paste0(gwas_table$file_name[i], "_processedgwas.tsv"))
write.table(snps_df2, file = outputfile, quote = F, row.names = F)

#Also write out a file with just the snps we are interested in for FGFP
outputfile_snps <- file.path(outputdir, paste0(gwas_table$file_name[i], "_snps.txt"))
write.table(snps_df2["fgfp_linker"], file = outputfile_snps, col.names=F, row.names=F, quote=F)
