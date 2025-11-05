#### We also want to look at some MiBioGen sumstats to see if we are observing the same pattern

#Lib
library("dplyr")
library("data.table")
library("RColorBrewer")
library("purrr")
library("TwoSampleMR")
library("ggplot2")

#Data path for proj dir
data.path <- "../../../data/PhenoPRS"

#Also load in path to scratch, proj_dir and fgfp data
source("../../parameters/base_dir.R")
source("../../parameters/scratch_path.txt")

#Also load in the SNP data
load(file.path(data.path, "phenosumstats/filtered_sumstats_list.RData"))

#Load in the candidate bugs
load(file.path(data.path, "data_out/data_ConsistentAssoc03.2.RData"))

#Load in our MBG trait SNPs
load(file.path(data.path, "data_out/BetaCor_05.3.RData"))

#Rename snps file with the new naming
#Specify a function for cleaning directory names
clean_directory_name <- function(x) {
  x %>%
    stringr::str_replace_all("[^[:alnum:]]+", "_") %>% 
    stringr::str_replace_all("_+", "_")   
}


names(snps_df) <- clean_directory_name(names(snps_df))


####################################
####### Load in LD SNPs ############
####################################

#Get the LD SNPs from our covariates of interest - lets start with r >= 0.1
ld_snps_out <- file.path(data.path, "phenosumstats/ldsnps4alphadiv.RData")

covars <- candidate_bugs$prs %>% unique()

if(!file.exists(ld_snps_out)){
  ld_snps <- list()
  for(covar in covars){
    ld_snps[[covar]] <- read.table(file.path(data.path, "phenosumstats", covar, "corsnps", paste0("corsnps_", covar, ".ld")), header=T) %>%
      filter(R2 >= 0.1)
  }
  
  #Save this is an .RData file so we can load locally
  save(ld_snps, file = ld_snps_out)
}else{load(ld_snps_out)}

ld_snps_long <- lapply(names(ld_snps), function(x) {
  cbind(ld_snps[[x]], x) %>% as.data.frame()
}) %>% do.call(rbind, .)


#########################################
# Make Manhattan plot for all FGFP traits
#########################################


#Remove urea as no longer in the final analysis
#ld_snps_long <- filter(ld_snps_long, x != "Serum_urea_levels")

#Focus on SNPs in high LD
ld_snps_long_highr2 <- filter(ld_snps_long, R2 > 0.8)

#Edit naming so nice for display
ld_snps_long_highr2 <- ld_snps_long_highr2 %>%
  mutate(covar_clean = case_when(
    x == "yengo_bmi" ~ "BMI",
    x == "Alanine_aminotransferase_levels" ~ "AAT",
    x == "Gamma_glutamyltransferase_levels_UKB_data_field_30730_" ~ "GGT",
    x == "Hip_circumference_UKB_data_field_49_" ~ "Hip circumference",
    x == "Standing_height" ~ "Standing height",
    x == "Red_blood_cell_count" ~ "RBC count",
    x == "HDL_cholesterol" ~ "HDL",
    x == "Serum_uric_acid_levels" ~ "Uric Acid",
    x == "triglycerides" ~ "Triglycerides",
    x == "Serum_urea_levels" ~ "Urea"
  ))


#Get vector of all our FGFP traits
fgfptraits_all <- unique(lm_out$bug)

#We need to write a function which reads in the FGFP sumstat of interest and then filters to p < 1e-4
#First make df to store results

load(file.path(data.path, "davesRNTsumstats", "fgfpfilt_sumstats.RData")) 

bug_sumstats <- out[names(out)[grepl("RNT", names(out))]] 


#Lets collapse into a single df
sigsnps_df <- lapply(names(bug_sumstats), function(x) {
  cbind(filter(bug_sumstats[[x]], frequentist_add_pvalue < 1e-4), bac = x)
}) %>% do.call(rbind, .)

allsnps_df <- lapply(names(bug_sumstats), function(x) {
  cbind(bug_sumstats[[x]], bac = x)
}) %>% do.call(rbind, .)

#Now filter to the SNPs in the MBG sumstats which are also suggestive
tofollow_df <- filter(sigsnps_df, rsid %in% ld_snps_long_highr2$SNP_B) %>% merge(ld_snps_long_highr2, by.x = "rsid", by.y = "SNP_B", all.x = T)

#Set our colours as the plot for consistency
traits <- tofollow_df$covar_clean %>% unique() 

#Remove blue colour as doesn't work with the blue Man plot
#colours <- brewer.pal(n = length(traits) + 1, name = "Spectral")[!brewer.pal(n = length(traits) + 1, name = "Spectral") %in% "#3288BD" ]

colours <- c("#E63946", "#2A9D8F", "#F4A261", "#9B5DE5", "#FF577F", "#4CE0D2", "#FFBE0B", "#6A4C93",  "#90BE6D", "#F9C74F", "#F3722C", "#F8961E", "#90BE6D", "#277DA1")

custom_colours <- setNames(colours, traits)


#Now make a df with the full names of our traits so we can add this to our plots
#full_names <- filter(lm_out_mbg, mbgbugs %in% unique(mbgtofollow_df$bac)) %>% select(mbgbugs, Phylum, Class, Order, Family, Genus) %>% unique()
#
##Remove those which we couldn't resolve to the same taxonomic level as FGFP i.e. duplicated
#full_names <- filter(full_names, !Genus %in% c("Sporobacter", "Coprococcus"))
#
#full_names$fullname <- apply(full_names[,2:6], 1, function(x){
#  x %>% paste(collapse = "_")
#})




#### Now make Man plot

#Lets set our x axis the same for every plot
combined_df <- allsnps_df %>% 
  select(chromosome, position) %>% 
  distinct()

# Compute cumulative positions
global_data_cum <- combined_df %>%
  group_by(chromosome) %>%
  summarise(max_bp = max(position)) %>%
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>%
  select(chromosome, bp_add, max_bp)

#Set an empty list for our manplots
manplot <- list()

matchedsnps_frommanplots <- list()

## Now we can loop round to make our manhattan plots for the traits of interest
for(trait in unique(tofollow_df$bac)){
  
  #Extract the full names for annotation
  #mt_fullname <- filter(full_names, mbgbugs == trait) %>% pull(fullname)
  
  snpsofinterest <- filter(tofollow_df, bac %in% trait)
  
  split_snps <- split(snpsofinterest, snpsofinterest$SNP_A)
  
  matchedsnps_frommanplots[[trait]] <- split_snps
  
  ################################
  ###### Visualisation of results
  ################################
  
  
  df <- filter(allsnps_df, bac %in% trait) %>% select(chromosome, frequentist_add_pvalue, rsid, position)
  
  names(df) <- c("CHR", "P", "SNP", "BP")
  
  gwas_data <- df %>%
    inner_join(global_data_cum, by = c("CHR" = "chromosome")) %>%
    mutate(bp_cum = BP + bp_add)
  
  axis_set <- global_data_cum %>%
    group_by(chromosome) %>%
    summarize(center = mean(bp_add + max_bp / 2))
  
  ylim <- gwas_data %>%
    filter(P == min(P)) %>%
    mutate(ylim = abs(floor(log10(P))) + 2) %>%
    pull(ylim)
  
  
  highlight_df <- snpsofinterest %>%
    mutate(R2 = round(R2, 2)) %>%
    merge(gwas_data, by.x = "rsid", by.y = "SNP", all.x = TRUE)
  
  
  #Add in the R2 annotation on the plot
  #Add in the R2
  annot_r2 <- highlight_df %>%
    group_by(SNP_A) %>%
    reframe(max_value = max(R2),
            bp_cum = median(bp_cum),
            y = max(-log10(P)),
            trait = unique(covar_clean),
            chr = unique(chromosome))
  
  annot_r2 <- annot_r2 %>%
    mutate(
      chr_group = ifelse(chr %in% c(9, 10), "9_10", as.character(chr)) # Group 9 and 10 together
    ) %>%
    arrange(bp_cum) %>% # Sort by x-axis positions
    group_by(chr_group) %>% # Group by chr_group
    mutate(
      y_offset = ifelse(
        chr_group == "9_10", 
        -log10(1e-2) + 0.4 * (row_number() - 1), # For 9 and 10 group
        -log10(1e-3) + 0.4 * (row_number() - 1)  # For other groups
      )
    )
  
  
  
  manplot[[trait]] <- ggplot(gwas_data, aes(x = bp_cum, y = -log10(P))) +
    # Genome-wide significance lines
    geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dashed", linewidth = 0.8) +
    geom_hline(yintercept = -log10(1e-5), color = "orange", linetype = "dashed", linewidth = 0.8) +
    
    # Manhattan dots
    geom_point(color = "#276FBF", alpha = 0.8, size = 1.2) +
    
    # Highlight specific SNPs with unified color mapping
    geom_point(data = highlight_df, aes(x = bp_cum, y = -log10(P), color = covar_clean), 
               size = 7, shape = 17) +  # Triangle shape
    
    # Optional annotations (uncomment if you want text labels)
    # geom_text(data = annot_r2, aes(x = bp_cum, y = y_offset, label = paste0("R² = ", max_value), color = trait),
    #           size = 3, angle = 45, vjust = -0.5, hjust = 0.5, show.legend = FALSE) +
    
    # Custom color mapping
    scale_color_manual(values = custom_colours) +
    
    # Axis labels and limits
    scale_x_continuous(label = axis_set$chromosome, breaks = axis_set$center) +
    scale_y_continuous(limits = c(-log10(1e-4), ylim), expand = c(0, 0)) +
    
    # Axis titles and plot title
    labs(
      x = "Chromosome",
      y = expression(-log[10](italic(p))),
      colour = "Covariate",
      title = paste0("Manhattan Plot of FGFP Trait: ", trait)
      #subtitle = paste0("Full Taxa Name: ", mt_fullname)
    ) +
    
    # Annotate significance thresholds
    annotate("text", 
             x = global_data_cum[[21,2]], 
             y = -log10(5e-8), 
             label = "p < 5e-8", 
             color = "red", 
             hjust = -0.1, 
             vjust = 0, 
             size = 3) +
    annotate("text", 
             x = global_data_cum[[21,2]], 
             y = -log10(1e-5), 
             label = "p < 1e-5", 
             color = "orange", 
             hjust = -0.1, 
             vjust = 0, 
             size = 3) +
    
    # Publication-style theme
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 60, size = 9, vjust = 1),
      axis.title = element_text(face = "bold", size = 13),
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  
}


##################################################################################################
##################################################################################################
#Now lets further look at the covariable SNPs proping up

#First we want to check where there are blocks of Covar in one region if they are picking up the same SNP
#lapply(matchedsnps_frommanplots[["family.Lachnospiraceae.id.1987"]], function(x) {
#  x$rsID}) %>% unlist() %>% table()
#Yes they are


#Now lets get a vector of all the group of SNP names
covarsnps_tofolow <- lapply(matchedsnps_frommanplots, function(x) {
  names(x) }) %>% stack() %>% rename(SNP = values, mbgbug = ind)

#We can also extract the block of SNPs relating to each SNP region of interest and the covar name
result_df <- map_dfr(matchedsnps_frommanplots, function(x) {
  map_dfr(x, function(y) {
    tibble(
      covar_name = unique(y$x),
      covar_clean = unique(y$covar_clean),
      rsID_block = paste(y$rsid, collapse = ",")
    )
  })
})


covarsnps_tofolow <- cbind(covarsnps_tofolow, result_df)




#Lets make a output list to hold our results
snpstofollow_results <- list()

#Change the factor variables to chr in the df
covarsnps_tofolow$mbgbug <- as.character(covarsnps_tofolow$mbgbug)

fullsspath <- file.path(data.path, "davesRNTsumstats")

#Loop round bugs
for(loop in 1:nrow(covarsnps_tofolow)){
  bug4loop <- covarsnps_tofolow[loop, "mbgbug"]
  i <- covarsnps_tofolow[loop, "covar_name"]
  
  snpstofollow_results[[loop]] <- list()
    
    #First get the singluar SNPs
    snps_vector <- covarsnps_tofolow[loop, "SNP"]
    
    #Then extract our list of SNPs
    snps_vector <- c(snps_vector,
                     filter(covarsnps_tofolow, mbgbug %in% bug4loop & covar_name == i) %>% pull(rsID_block) %>% strsplit(",") %>% unlist() %>% unique())
    
    #Now we want the SNPs from bug GWAS
    snpstofollow_results[[loop]][["bac"]] <- 
      tofollow_df %>% 
      filter(rsid %in% snps_vector & x == i)
    
    #And we also want the SNP results from the covar GWAS 
    snpstofollow_results[[loop]][["covar"]] <- 
      filter(snps_df[[i]], rsid %in% covarsnps_tofolow[loop, "SNP"])
    
}

covarsnps_tofollow <- cbind(covarsnps_tofolow,
                            map_dfr(snpstofollow_results, function(x) {
                                x$covar
                              })
)


#Lets have a look at these covar SNPs to see if they are strong associations
#As these will just be the one SNP, we can make this into a table


#Now we want to look at our bug SNPs and add them to the dataframe
#First we need to find these SNPs in full sumstats to make sure gwas data is harmonised, so no confusion with directionality



snpsearch<- NULL

for(i in 1:nrow(covarsnps_tofollow)){
  snpsearch <- rbind(snpsearch,
                     data.table::fread(file = file.path(fullsspath,
                                                        paste0(covarsnps_tofollow[[i, "mbgbug"]], "_allchr.txt.gz"))) %>% 
                       filter(rsid %in% covarsnps_tofollow[[i, "rsid"]])
  )
}

colnames(snpsearch) <- paste0(colnames(snpsearch), "_bug")

#Lets bind this output to our covarsnps_tofollow df
covarsnps_tofollow <- cbind(covarsnps_tofollow, snpsearch)

# Now we can harmonise our two datasources
# First check if any of these are pallindromic and will require further data for resolution
covarsnps_tofollow <- covarsnps_tofollow %>%
  mutate(
    # Identify palindromic SNPs
    is_palindromic = (alleleA_bug == "A" & alleleB_bug == "T") |
      (alleleA_bug == "T" & alleleB_bug == "A") |
      (alleleA_bug == "C" & alleleB_bug == "G") |
      (alleleA_bug == "G" & alleleB_bug == "C"))

stopifnot(sum(covarsnps_tofollow$is_palindromic) == 0)

#We have no palindromic SNPs

# Now check that we have the same alleles for the same SNPs
stopifnot(apply(covarsnps_tofollow, 1, function(x) {
  !(x["eff.allele"] %in% x[c("alleleA_bug", "alleleB_bug")])
}) %>% sum() ==0)

#Now we can flip if effect allele is wrong
covarsnps_tofollow_flipped <- covarsnps_tofollow %>% mutate(
  toflip = ifelse(alleleB_bug == eff.allele, F, T)
) %>%
  #Now we can flip alleles and beta dir if Y
  mutate(
    eff.allele = alleleB_bug,
    beta = ifelse(toflip, -beta, beta),
    EAF = ifelse(toflip, 1 - EAF, EAF)
  ) %>%
  #Lets also remove the Oth_al column as most of this data is missing and hasn't been corrected for the flipping
  select(-Oth_al)

#Check our alleles are all the same now
stopifnot(sum(covarsnps_tofollow_flipped$alleleB_bug != covarsnps_tofollow_flipped$alleleB_bug) == 0)

#Change colnames so no confusion between the two datasources
colnames(covarsnps_tofollow_flipped)[6:11] <- 
  paste0(colnames(covarsnps_tofollow_flipped)[6:11], "_covar")

#Now we can merge with the SNP with lowest pval from bug sumstats so we can visualise our results
finaldf_snpstofollow <- do.call(rbind, lapply(seq_along(snpstofollow_results), function(i) {
  x <- snpstofollow_results[[i]]
  
    filtered <- x$bac %>% filter(frequentist_add_pvalue == min(frequentist_add_pvalue)) %>%
      rename(covar_name = x)
    return(filtered)
  
  
}))  %>%
  cbind(covarsnps_tofollow_flipped)


finaldf_snpstofollow <- as.data.frame(finaldf_snpstofollow)

clean_finaldf <- finaldf_snpstofollow[,which(!duplicated(names(finaldf_snpstofollow)))] %>%
  select(
    #First select the lowest pval from the region
    mbgbug, rsid, chromosome, position, frequentist_add_pvalue, covar_clean, covar_name, 
    #Now select the covar SNP and R2 between the two SNPs
    SNP, R2,
    #Now the results of the bug GWAS for this SNP
    alleleB_bug, frequentist_add_beta_1_bug, frequentist_add_se_1_bug,
    #And finally the covar SNP where it was discovered
    eff.allele_covar, beta_covar, SE_covar, pval_covar
  )

snpstofollow_dirs <- list(snpstofollow_dirs = finaldf_snpstofollow, clean_snpstofollow = clean_finaldf)

#Save this output so can be called upon in the Markdown
save(manplot, matchedsnps_frommanplots, snpstofollow_dirs,
     file = file.path(data.path, "data_out/ManPlotsFGFP_05.6.RData"))




###############
# We also want to make a full manhattan plot for G_unclassified_P_Firmicutes_RNT as this is important step in our analyses
#Lets read in the full summary statistics for the FGFP Microbiome result
fullsspath <- file.path(data.path, "davesRNTsumstats")
fullss <- fread(file.path(fullsspath, "G_unclassified_P_Firmicutes_RNT_allchr.txt.gz"))

highlight_snps <- lapply(matchedsnps_frommanplots[["G_unclassified_P_Firmicutes_RNT"]], function(x) x$rsid) %>% unlist(use.names = F)

manplot <- manhattan(fullss,
                     chr = "chromosome",
                     bp = "position",
                     snp = "rsid",
                     p = "frequentist_add_pvalue", 
                     highlight = highlight_snps) 


manplot
