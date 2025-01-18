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
load(file.path(data.path, "data_out/BetaCorMBG_05.4.RData"))

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
# Make Manhattan plot for all MBG traits
#########################################


#Remove urea as no longer in the final analysis
ld_snps_long <- filter(ld_snps_long, x != "Serum_urea_levels")

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
    x == "triglycerides" ~ "Triglycerides"
  ))


#First check if any alpha_div SNPs show up as these are in a difference datasource
alphadivsnps <- read.csv(file = file.path(data.path, "mibiogen_ss/mbg_alphadivSNPs.csv"))

filter(ld_snps_long_highr2, SNP_B %in% alphadivsnps$rsID) 
#Nothing in alpha diversity

#Get vector of all our MBG traits
mbgtraits_all <- unique(lm_out_mbg$mbgbugs)

#Load in the MGBe4_snps
load(file.path(data.path, "mibiogen_ss/MBGe4_snps.RData"))

mbg_snps <- filter(MBGe4_snps, bac %in% mbgtraits_all)

#Now filter to the SNPs in the MBG sumstats which are also suggestive
mbgtofollow_df <- filter(mbg_snps, rsID %in% ld_snps_long_highr2$SNP_B) %>% merge(ld_snps_long_highr2, by.x = "rsID", by.y = "SNP_B", all.x = T)

#Set our colours as the plot for consistency
traits <- mbgtofollow_df$covar_clean %>% unique() 

#Remove blue colour as doesn't work with the blue Man plot
#colours <- brewer.pal(n = length(traits) + 1, name = "Spectral")[!brewer.pal(n = length(traits) + 1, name = "Spectral") %in% "#3288BD" ]

colours <- c("#E63946", "#2A9D8F", "#F4A261", "#9B5DE5", "#FF577F", "#4CE0D2", "#FFBE0B", "#6A4C93")

custom_colours <- setNames(colours, traits)


#Now make a df with the full names of our traits so we can add this to our plots
full_names <- filter(lm_out_mbg, mbgbugs %in% unique(mbgtofollow_df$bac)) %>% select(mbgbugs, Phylum, Class, Order, Family, Genus) %>% unique()

#Remove those which we couldn't resolve to the same taxonomic level as FGFP i.e. duplicated
full_names <- filter(full_names, !Genus %in% c("Sporobacter", "Coprococcus"))

full_names$fullname <- apply(full_names[,2:6], 1, function(x){
  x %>% paste(collapse = "_")
})




#### Now make Man plot

#Lets set our x axis the same for every plot
combined_df <- MBGe4_snps %>% 
  select(chr, bp) %>% 
  distinct()

# Compute cumulative positions
global_data_cum <- combined_df %>%
  group_by(chr) %>%
  summarise(max_bp = max(bp)) %>%
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>%
  select(chr, bp_add, max_bp)

#Set an empty list for our manplots
manplot <- list()

matchedsnps_frommanplots <- list()

## Now we can loop round to make our manhattan plots for the traits of interest
for(trait in unique(mbgtofollow_df$bac)){
  
  #Extract the full names for annotation
  mt_fullname <- filter(full_names, mbgbugs == trait) %>% pull(fullname)
  
  snpsofinterest <- filter(mbgtofollow_df, bac %in% trait)
  
  split_snps <- split(snpsofinterest, snpsofinterest$SNP_A)
  
  matchedsnps_frommanplots[[trait]] <- split_snps
  
  ################################
  ###### Visualisation of results
  ################################
  
  
  df <- filter(MBGe4_snps, bac %in% trait) %>% select(chr, P.weightedSumZ, rsID, bp)
  
  names(df) <- c("CHR", "P", "SNP", "BP")
  
  gwas_data <- df %>%
    inner_join(global_data_cum, by = c("CHR" = "chr")) %>%
    mutate(bp_cum = BP + bp_add)
  
  axis_set <- global_data_cum %>%
    group_by(chr) %>%
    summarize(center = mean(bp_add + max_bp / 2))
  
  ylim <- gwas_data %>%
    filter(P == min(P)) %>%
    mutate(ylim = abs(floor(log10(P))) + 2) %>%
    pull(ylim)
  
  
  highlight_df <- snpsofinterest %>%
    mutate(R2 = round(R2, 2)) %>%
    merge(gwas_data, by.x = "rsID", by.y = "SNP", all.x = TRUE)
  
  
  #Add in the R2 annotation on the plot
  #Add in the R2
  annot_r2 <- highlight_df %>%
    group_by(SNP_A) %>%
    reframe(max_value = max(R2),
            bp_cum = median(bp_cum),
            y = max(-log10(P)),
            trait = unique(covar_clean),
            chr = unique(chr))
  
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
  
  
  
  manplot[[trait]] <- ggplot(gwas_data, aes(
    x = bp_cum, 
    y = -log10(P)
  )) +
    geom_hline(
      yintercept = -log10(5e-8), 
      color = "red", 
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(1e-5), 
      color = "red", 
      linetype = "dashed"
    ) +
    # Add shaded area from 0 to -log10(1e-4) across the x-axis
    geom_rect(
      xmin = min(gwas_data$bp_cum), xmax = max(gwas_data$bp_cum), 
      ymin = 0, ymax = -log10(1e-4), 
      fill = "gray", alpha = 0.3
    ) +
    # Set points to a consistent grey color
    geom_point(color = "#276FBF", alpha = 0.75, size = 1.5) +
    
    # Highlight specific SNPs with unified color mapping
    geom_point(data = highlight_df, aes(
      x = bp_cum, 
      y = -log10(P.weightedSumZ), 
      color = covar_clean,
    ), size = 5, shape = 17) + # Triangles
    
    # # Add annotations with unified color mapping
    # geom_text(
    #   data = annot_r2,
    #   aes(
    #     x = bp_cum, 
    #     y = y_offset, 
    #     label = paste0("R2 = ", max_value),
    #     color = trait  
    #   ),
    #   vjust = -0.5,     
    #   hjust = 0.5,      
    #   size = 3,
    #   angle = 45,
    #   show.legend = F
    # ) +
    
    # Define a single color scale for both triangles and annotations
    scale_color_manual(values = custom_colours) +
    
    # Add annotation for p < 5e-8 line
    annotate("text", 
             x = global_data_cum[[21,2]], 
             y = -log10(5e-8), 
             label = "p < 5e-8", 
             color = "red", 
             hjust = -0.1, 
             vjust = 0) +
    # Add annotation for p < 1e-5 line
    annotate("text", 
             x = global_data_cum[[21,2]], 
             y = -log10(1e-5), 
             label = "p < 1e-5", 
             color = "red", 
             hjust = -0.1, 
             vjust = 0) +
    
    # Formatting for the axes and theme
    scale_x_continuous(
      label = axis_set$chr, 
      breaks = axis_set$center
    ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +
    labs(
      x = "CHR", 
      y = "-log10(p)",
      colour = "Covariate",
      title = paste0("Manhattan Plot of MiBioGen Bug: ", trait),
      subtitle = paste0("Full Taxa Name: ", mt_fullname)
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 60, size = 8, vjust = 2), 
      plot.title = element_text(size = 20, face = "bold")
    )
  
}


##################################################################################################
##################################################################################################
#Now lets further look at the covariable SNPs proping up

#First we want to check where there are blocks of Covar in one region if they are picking up the same SNP
lapply(matchedsnps_frommanplots[["family.Lachnospiraceae.id.1987"]], function(x) {
  x$rsID}) %>% unlist() %>% table()
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
      rsID_block = paste(y$rsID, collapse = ",")
    )
  })
})


covarsnps_tofolow <- cbind(covarsnps_tofolow, result_df)




#Lets make a output list to hold our results
snpstofollow_results <- list()

#Change the factor variables to chr in the df
covarsnps_tofolow$mbgbug <- as.character(covarsnps_tofolow$mbgbug)



#Loop round bugs
for(bug4loop in unique(covarsnps_tofolow$mbgbug)){
#And loop round the covars
for(i in filter(covarsnps_tofolow, mbgbug %in% bug4loop)[["covar_name"]]){

#First get the singluar SNPs
snps_vector <- filter(covarsnps_tofolow, mbgbug %in% bug4loop & covar_name == i) %>%
  pull(SNP)

#Then extract our list of SNPs
snps_vector <- c(snps_vector,
                 filter(covarsnps_tofolow, mbgbug %in% bug4loop & covar_name == i) %>% pull(rsID_block) %>% strsplit(",") %>% unlist() %>% unique())

#Now we want the SNPs from bug GWAS
snpstofollow_results[[bug4loop]][[i]][["bac"]] <- 
  filter(MBGe4_snps, bac %in% bug4loop, rsID %in% snps_vector)

#And we also want the SNP results from the covar GWAS 
snpstofollow_results[[bug4loop]][[i]][["covar"]] <- 
  filter(snps_df[[i]], rsid %in% snps_vector)

}
}

covarsnps_tofollow <- cbind(covarsnps_tofolow,
                            map_dfr(snpstofollow_results, function(x) {
                              map_dfr(x,function(y){
                                y$covar
                              })
                            })
)


#Lets have a look at these covar SNPs to see if they are strong associations
#As these will just be the one SNP, we can make this into a table


#Now we want to look at our bug SNPs and add them to the dataframe
#First we need to find these SNPs in full sumstats to make sure gwas data is harmonised, so no confusion with directionality

fullsspath <- file.path(data.path, "mibiogen_ss")

snpsearch<- NULL

for(i in 1:nrow(covarsnps_tofollow)){
  snpsearch <- rbind(snpsearch,
                     data.table::fread(file = file.path(fullsspath,
                                                  paste0(covarsnps_tofollow[[i, "mbgbug"]], ".summary.txt.gz"))) %>% 
    filter(rsID %in% covarsnps_tofollow[[i, "rsid"]])
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
    is_palindromic = (ref.allele_bug == "A" & eff.allele_bug == "T") |
      (ref.allele_bug == "T" & eff.allele_bug == "A") |
      (ref.allele_bug == "C" & eff.allele_bug == "G") |
      (ref.allele_bug == "G" & eff.allele_bug == "C"))

stopifnot(sum(covarsnps_tofollow$is_palindromic) == 0)

#We have no palindromic SNPs

# Now check that we have the same alleles for the same SNPs
stopifnot(apply(covarsnps_tofollow, 1, function(x) {
  !(x["eff.allele"] %in% x[c("eff.allele_bug", "ref.allele_bug")])
}) %>% sum() ==0)

#Now we can flip if effect allele is wrong
covarsnps_tofollow_flipped <- covarsnps_tofollow %>% mutate(
  toflip = ifelse(eff.allele_bug == eff.allele, F, T)
) %>%
  #Now we can flip alleles and beta dir if Y
  mutate(
    eff.allele = eff.allele_bug,
    beta = ifelse(toflip, -beta, beta),
    EAF = ifelse(toflip, 1 - EAF, EAF)
  ) %>%
  #Lets also remove the Oth_al column as most of this data is missing and hasn't been corrected for the flipping
  select(-Oth_al)

#Check our alleles are all the same now
stopifnot(sum(covarsnps_tofollow_flipped$eff.allele_bug != covarsnps_tofollow_flipped$eff.allele) == 0)

#Change colnames so no confusion between the two datasources
colnames(covarsnps_tofollow_flipped)[6:11] <- 
  paste0(colnames(covarsnps_tofollow_flipped)[6:11], "_covar")

#Now we can merge with the SNP with lowest pval from bug sumstats so we can visualise our results
finaldf_snpstofollow <- map_dfr(snpstofollow_results, function(x) {
  map_dfr(x, function(y) {
    filter(y$bac, P.weightedSumZ == min(P.weightedSumZ)) 
  }, .id = "covar_name")
}) %>% cbind(covarsnps_tofollow_flipped)

#Find R2 of the SNP
finaldf_snpstofollow$R2 <- 
  apply(finaldf_snpstofollow, 1, function(df) {
  filt <- dplyr::filter(ld_snps_long_highr2, x == df["covar_name"] & SNP_B == df["rsID"] & SNP_A == df["rsid_covar"]) %>% pull(R2)
  if(length(filt) == 0){
    return(NA)}else{return(filt)}
}) %>% unlist()


clean_finaldf <- finaldf_snpstofollow[,!duplicated(names(finaldf_snpstofollow))] %>%
  select(
    #First select the lowest pval from the region
    bac, rsID, chr, bp, P.weightedSumZ, N, covar_clean, covar_name, 
    #Now select the covar SNP and R2 between the two SNPs
    SNP, R2,
    #Now the results of the bug GWAS for this SNP
    eff.allele_bug, beta_bug, SE_bug,
    #And finally the covar SNP where it was discovered
    eff.allele_covar, beta_covar, SE_covar, pval_covar
    )

snpstofollow_dirs <- list(snpstofollow_dirs = finaldf_snpstofollow, clean_snpstofollow = clean_finaldf)

#Save this output so can be called upon in the Markdown
save(manplot, matchedsnps_frommanplots, snpstofollow_dirs,
     file = file.path(data.path, "data_out/ManPlotsMBG_05.5.RData"))

