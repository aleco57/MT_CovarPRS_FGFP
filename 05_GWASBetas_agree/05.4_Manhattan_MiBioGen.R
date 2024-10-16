#### We also want to look at some MiBioGen sumstats to see if we are observing the same pattern

#### MiBioGen
################
#Lets test with the Roseburia bug
#Read in the sumstat
sumstat <- fread(file.path(data.path, "mibiogen_ss/genus.Roseburia.id.2012.summary.txt.gz"))

#Which SNPs are we intersted in? 
filter(candidate_bugs, mt == "G_Roseburia_RNTRes")

snps <- c(snps_df[["Alanine aminotransferase levels"]]$rsid, snps_df[["yengo_bmi"]]$rsid)

filt_ss <- filter(sumstat, rsID %in% snps)

#Do we have signal here
filter(filt_ss, P.weightedSumZ < 1e-2)

#No signal here

