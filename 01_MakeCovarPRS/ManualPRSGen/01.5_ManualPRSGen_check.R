## Set wd to PhenoPRS
#Set up
library(dplyr)

#Data path for proj dir
data.path <- "../../../data/PhenoPRS"

#Set trait variable - this can also be changed to read in from bash script
trait <- "creatinine3"
trait_path <- file.path(data.path, "prs_out", trait)

#read in gen and sample file
gen <- read.table(file.path(trait_path, paste0(trait, ".gen")), as.is =T)
sample <- read.table(file.path(trait_path, paste0(trait, ".sample")), as.is =T, skip = 2)
  
#Read in beta GWAS file
gwas <- read.table(file.path(trait_path, paste0(trait, "_processedgwas.tsv")), as.is =T, sep = " ", header = T)

#First we can make our function to calculate dosages from the file

#Num_fixed cols in the non dosage part at beginning of .gen file
calc_dosage <- function(gen_data, num_fixed_cols) {
  num_individuals <- (ncol(gen_data) - num_fixed_cols) / 3
  dosage_matrix <- matrix(NA, nrow = nrow(gen_data), ncol = num_individuals)
  
  for (i in 1:num_individuals) {
    pAA_col <- num_fixed_cols + (i - 1) * 3 + 1
    pAB_col <- pAA_col + 1
    pBB_col <- pAB_col + 1
    
    # Calculate dosage of allele B
    dosage <- 0 * gen_data[[pAA_col]] + 1 * gen_data[[pAB_col]] + 2 * gen_data[[pBB_col]]
    
    # Store the dosage in the matrix
    dosage_matrix[, i] <- dosage
  }
  
  #Return with the other parts of the df
  return(cbind(gen_data[,1:num_fixed_cols],dosage_matrix))
}

dosages <- calc_dosage(gen,6)

#Merge with beta file
dosage_beta <- merge(gwas, dosages, by.x = "fgfp_linker", by.y = "V3")

#Does our eff.allele correspond to our dosage?
identical(dosage_beta$eff.allele, dosage_beta$V6)
#Yes it does

#We can make our score in two ways, first with the beta from file and second with only + beta
#Lets now first generate out score using beta values which are there
prs1 <- colSums(dosage_beta[, 13:ncol(dosage_beta)] * dosage_beta$beta) / 2*nrow(dosage_beta) #Here we are dividing by number of alleles to standardise
#Now we can change so that everything is on a positive scale, meaning easier interpretation
#Here we chanfe the dosage 
dosage_beta_pos <- dosage_beta
dosage_beta_pos[dosage_beta$beta < 0, 13:ncol(dosage_beta)] <- 2 - (dosage_beta[dosage_beta$beta < 0, 13:ncol(dosage_beta)])
dosage_beta_pos$beta <- abs(dosage_beta_pos$beta)
prs2 <- colSums(dosage_beta_pos[, 13:ncol(dosage_beta_pos)] * dosage_beta_pos$beta) / 2*nrow(dosage_beta_pos)

#Now we can generate our PRS file and scale them

prs <- data.frame(ID = sample$V1, prs1 = prs1, prs_posbeta = prs2) #Here we are dividing by number of alleles to standardise

#NOTE!!! When we scale and centre the two PRSs generated become identical... obviously

#Save our PRS
write.table(prs, file = file.path(data.path, "prs_out", trait, paste0(trait, "_manualPRS_gen.txt")), quote =F, row.names = F)

#Load in phenotypic data to see how our PRSs compare
load(file.path(data.path, "fgfpdata4prs.RData"))

prs_trait <- merge(prs,
                   data4prs[["pheno_covariate_prs"]][,c("IID", "creatinine3.Pt_5e.08","Creatinine_mgdL", "e.GFR")] %>% distinct(),
                   by.x = "ID",
                   by.y = "IID",
                   all.x = T)

apply(prs_trait[,2:4], 2, scale)   
#We can see when scaling the result for creatinine generated from PRSice-2 is slightly out
## NOTE! This could potnentially be due to handeling of missing data via MAF impute !! or if not doing dosage calculation for PRS

#Now lets plot scaled and non-scaled plots of our creatinine PRS

#First we can look at the scatter plots
ggplot(prs_trait, aes(x = prs_posbeta, y = Creatinine_mgdL)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", color = "red") +
  labs(title = "Scatter Plot with Regression Line",
       x = "Manual PRS of Creatinine", y = "Creatinine_mgdL")


#We can load in GFR from the data and see how this associates
#First edit the phenotype compile data to add eGFR
ggplot(prs_trait, aes(x = prs_posbeta, y = e.GFR)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", color = "red") +
  labs(title = "Scatter Plot with Regression Line",
       x = "Manual PRS of Creatinine", y = "e.GFR")
