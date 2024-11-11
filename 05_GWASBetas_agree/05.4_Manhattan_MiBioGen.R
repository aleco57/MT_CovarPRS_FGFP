#### We also want to look at some MiBioGen sumstats to see if we are observing the same pattern

#Lib
library("dplyr")
library("data.table")

#Data path for proj dir
data.path <- "../../../data/PhenoPRS"

#Also load in path to scratch, proj_dir and fgfp data
source("../../parameters/base_dir.R")
source("../../parameters/scratch_path.txt")

#Also load in the SNP data
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

#Are any of these SNPs showing up in the alpha diversity sumstats?
ld_snps_long <- lapply(names(ld_snps), function(x) {
  cbind(ld_snps[[x]], x) %>% as.data.frame()
}) %>% do.call(rbind, .)


###############################
# Individual bugs from MiBioGen
###############################

#Lets load in all the suggestive SNPs published by MiBioGen
load(file.path(data.path, "mibiogen_ss/MBGe4_snps.RData"))

MBGe4_snps <- read.table(file.path(data.path, "mibiogen_ss/MBG.allHits.p1e4.txt"), header=T, as.is = T)

#Which bugs can we look at?
unique(candidate_bugs$mt)

#filter(ld_snps_long, SNP_B %in% filter(MBGe4_snps, bac == "genus.Roseburia.id.2012" & rsID %in% snps4ros)[["rsID"]])

fgfpbugs <- c("G_Roseburia_RNTRes", 
              "F_Porphyromonadaceae_RNTRes",
              "G_Barnesiella_RNTRes",
              "F_Erysipelotrichaceae_RNTRes",
              "G_Ruminococcus_RNTRes",
              "G_Ruminococcus_RNTRes",
              "G_unclassified_F_Erysipelotrichaceae_RNTRes",
              "G_unclassified_F_Ruminococcaceae_RNTRes",
              "G_unclassified_P_Firmicutes_RNTRes",
              "F_Ruminococcaceae_RNTRes",
              "G_Oscillibacter_RNTRes",
              "G_Bacteroides_RNTRes",
              "G_Bifidobacterium_RNTRes")

mbgbugs <- c("genus.Roseburia.id.2012", 
             "family.Porphyromonadaceae.id.943",
             "genus.Barnesiella.id.944",
             "family.Erysipelotrichaceae.id.2149",
             "genus.Ruminococcus1.id.11373",
             "genus.Ruminococcus2.id.11374",
             "family.Erysipelotrichaceae.id.2149",
             "family.Ruminococcaceae.id.2050",
             "phylum.Firmicutes.id.1672",
             "family.Ruminococcaceae.id.2050",
             "genus.Oscillibacter.id.2063",
             "genus.Bacteroides.id.918",
             "genus.Bifidobacterium.id.436")

matched_bugs <- data.frame(fgfpbugs = fgfpbugs, mbgbugs = mbgbugs)


#For each match we will search for those covariate SNPs (and LDs) in the bug sumstats

results <- list()

for(i in 1:nrow(matched_bugs)){
  fgfpbug <- matched_bugs$fgfpbugs[i]
  mbgbug <- matched_bugs$mbgbugs[i]
  
  #Filter to the MBG associations for this specific bug
  mbg_snps <- filter(MBGe4_snps, bac == mbgbug)
  
  #Which covariates was this SNP associated with?
  covars <- filter(candidate_bugs, mt == fgfpbug)[["prs"]]
  snps <- filter(ld_snps_long, x %in% covars) %>% select(SNP_A, SNP_B) %>% unlist() %>% unique()
  
  #Now filter to the SNPs in the MBG sumstats
  mbg_snps <- filter(mbg_snps, rsID %in% snps)
  
  #Also want to merge the r2 values
  snpsofinterest <- merge(mbg_snps, 
                          filter(ld_snps_long, x %in% covars & SNP_B %in% mbg_snps$rsID) %>% select(SNP_A, SNP_B, R2, x),
                          by.x = "rsID", by.y = "SNP_B", all = T)
  
  results[[fgfpbug]] <- snpsofinterest
}

results_long <- lapply(names(results), function(x) {
  if(nrow(results[[x]]) > 0) {cbind(results[[x]], x) %>% as.data.frame()}
}) %>% do.call(rbind, .)
  
names(results_long)[15:16] <- c("covar", "fgfpbug")  


#Add in covaraite names which are better for visualisation / plotting
mutate(results_long, covar_clean = case_when(
  covar == "yengo_bmi" ~ "BMI",
  covar == "Alanine_aminotransferase_levels" ~ "AAT",
  covar == "Gamma_glutamyltransferase_levels_UKB_data_field_30730_" ~ "GGT",
  covar == "Hip_circumference_UKB_data_field_49_" ~ "Hip circumference",
  covar == "Standing_height" ~ "Standing height",
  covar == "Red_blood_cell_count" ~ "RBC count",
  covar == "Serum_urea_levels"  ~ "Serum Urea",
)) -> results_long


################################
###### Visualisation of results
################################

#Manhattan plots
manplots <- list()

#Make vector of those bugs with interesting plots
traits4plots <- c("genus.Roseburia.id.2012", 
                  "family.Porphyromonadaceae.id.943",
                  "family.Erysipelotrichaceae.id.2149",
                  "phylum.Firmicutes.id.1672",
                  "family.Ruminococcaceae.id.2050",
                  "genus.Bifidobacterium.id.436")


#Now we can assign a colour to each of the covariates on the plot, so the any one covar will be the same on diff plots
traits <- filter(results_long, bac %in% traits4plots)[["covar_clean"]] %>% unique()
custom_colours <- setNames(c("orange", "green", "purple", "maroon", "brown"), traits)

#Lets start with Rosaburia
for(mbgbug in traits4plots){

df <- filter(MBGe4_snps, bac == mbgbug) %>% select(chr, P.weightedSumZ, rsID, bp)

names(df) <- c("CHR", "P", "SNP", "BP")

#Make x axis for manhattan plot 
data_cum <- df %>%
  group_by(CHR) %>%
  summarise(max_bp = max(BP)) %>%
  mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>%
  select(CHR, bp_add)

gwas_data <- df %>%
  inner_join(data_cum, by = "CHR") %>%
  mutate(bp_cum = BP + bp_add)

axis_set <- gwas_data %>%
  group_by(CHR) %>%
  summarize(center = mean(bp_cum))

ylim <- gwas_data %>%
  filter(P == min(P)) %>%
  mutate(ylim = abs(floor(log10(P))) + 2) %>%
  pull(ylim)


highlight_df <- filter(results_long, bac == mbgbug) %>%
  mutate(R2 = round(R2, 2)) %>% merge(gwas_data, by.x = "rsID", by.y = "SNP", all.x = T)

#Add in the R2 annotation on the plot
#Add in the R2
annot_r2 <- highlight_df %>%
  group_by(SNP_A) %>%
  reframe(max_value = max(R2),
          bp_cum = median(bp_cum),
          y = max(-log10(P)),
          trait = unique(covar_clean))


manplots[[mbgbug]] <- ggplot(gwas_data, aes(
  x = bp_cum, 
  y = -log10(P)
)) +
  # Add shaded area from 0 to -log10(1e-4) across the x-axis
  geom_rect(
    xmin = min(gwas_data$bp_cum), xmax = max(gwas_data$bp_cum), 
    ymin = 0, ymax = -log10(1e-4), 
    fill = "gray", alpha = 0.3
  ) +
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
  # Set points to grey color
  geom_point(color = "#276FBF", alpha = 0.75, size = 1.5) +  # Grey points
  
  # Highlight specific SNPs from the new dataframe (colored by the 'x' column)
  geom_point(data = highlight_df, aes(
    x = bp_cum, y = -log10(P.weightedSumZ), color = covar_clean
  ), size = 3, shape = 17) +  # Changed shape to 17 (triangle) for visibility
  
  # Define your desired colors for specific values of 'x'
  scale_color_manual(values = custom_colours) + 
  # Add annotations using geom_text or geom_label
  geom_text(
    data = annot_r2,  # Use the annot_r2 dataframe for annotations
    aes(
      x = bp_cum,     # X coordinate from bp_cum
      y = 3,          # Y coordinate from y
      label = paste0(trait, " SNP: ", SNP_A, "\nSNP Region R2 = ", max_value)  # Annotation text
  ),
  vjust = -0.5,     # Adjust vertical position
  hjust = 0.5,      # Adjust horizontal position
  size = 3,
  angle = 45
  ) +
  # Add annotation for p < 5e-8 line
  annotate("text", 
           x = max(gwas_data$bp_cum),  # Place on the right side of the plot
           y = -log10(5e-8), 
           label = "p < 5e-8", 
           color = "red", 
           hjust = -0.1,  # Adjust for slight spacing from edge
           vjust = 0) +
  # Add annotation for p < 1e-5 line
  annotate("text", 
           x = max(gwas_data$bp_cum), 
           y = -log10(1e-5), 
           label = "p < 1e-5", 
           color = "red", 
           hjust = -0.1, 
           vjust = 0)  +
  
  # Formatting for the axes and theme
  scale_x_continuous(
    label = axis_set$CHR, 
    breaks = axis_set$center
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +
  
  labs(
    x = "CHR", 
    y = "-log10(p)",
  colour = "Covariate",
  title = paste0("Manhattan Plot of suggestive SNPs from MiBioGen Bug: ", mbgbug)) +
  theme_minimal() +
  theme(
    legend.position = "right",  # Position legend to the right
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5), 
    plot.title = element_text(size = 20, face = "bold")
  )

}


save(manplots, matched_bugs, file = file.path(data.path, "data_out/mibiogen_manplots.RData"))




#Now we can add our PRS / Observational forest plot above the manhattan plot for clarity of our results
#First reshape the data
# I think we should add in the beta from regression of sumstat + genetic cor where possible
obsest_prsest_long <- candidate_bugs %>%
  pivot_longer(cols = c(estimate_obs, estimate_prs), 
               names_to = "type", 
               values_to = "estimate") %>%
  pivot_longer(cols = c(std.error_obs, std.error_prs), 
               names_to = "error_type", 
               values_to = "std.error") %>%
  pivot_longer(cols = c(p.value_obs, p.value_prs), 
               names_to = "p_value_type", 
               values_to = "p.value") %>%
  filter((type == "estimate_obs" & error_type == "std.error_obs" & p_value_type == "p.value_obs") |
           (type == "estimate_prs" & error_type == "std.error_prs" & p_value_type == "p.value_prs")) %>%
  mutate(type = recode(type, estimate_obs = "Observational", estimate_prs = "PRS"))


#Now plot our traits of interest
#Forest plot for the observational and PRS estimate
ggforestplot::forestplot(
  df = filter(obsest_prsest_long, group != "Diversity"),
  name = mt,
  estimate = estimate,
  se = std.error,
  pvalue = p.value,
  colour = type,
  xlab = "SD increase in Rank normal transformed Microbial trait per 1-SD increase in polygenic loading of the covariate / the measure covariate",
  title = "Regression coefficents from Observational analyses of Covariate regressed on RNT microbial Traits in FGFP") +
  ggforce::facet_col(
    facets = ~covar_name,
    scales = "free_y"
  )




#Dir consistent? between covar GWAS and bug GWAS, 
##### !!!!we need MAF for this which is missing from mbg


####################################
##### alpha diversity from MiBioGen
####################################

#First read in the SNPs from supplementary MiBioGen table
alphadiv_snps <- read.csv(file.path(data.path, "mbg_alphadivSNPs.csv"))

#Now make a vector of covars that we are interested in from previous results
covars_alphadiv <- c("yengo_bmi","Alanine_aminotransferase_levels","triglycerides","Hip_circumference_UKB_data_field_49_","Serum_uric_acid_levels", "HDL_cholesterol")

#Make lapply to loop round all covars and get snps in one long df
snps_df_long <- lapply(covars, function(x) {
  cbind(snps_df[[x]]$rsid, x) %>% as.data.frame()
}) %>% do.call(rbind, .)

names(snps_df_long) <- c("rsid", "covar")

alphadiv_snps$rsID %in% snps_df_long$rsid %>% sum()

#No SNPs are in the initial, now lets look at SNPs also in LD with our instrument SNPs


correlatedsnps <- alphadiv_snps[alphadiv_snps$rsID %in% ld_snps_long$SNP_B %>% which(), ]

filter(ld_snps_long, SNP_B %in% correlatedsnps$rsID)

#We see an AAT correlated SNP in the Simpson GWAS results

#How many SNPs are suggesting in the Simpson results? 
filter(alphadiv_snps, ProbeName == "simpson")


#Lets look at the manhattan plot for this SNP
df <- filter(alphadiv_snps, ProbeName == "simpson") %>% select(chr, P, rsID, bp)
