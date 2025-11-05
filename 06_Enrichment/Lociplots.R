#### This will be an R script to make better Manhattan visualisation plots for the CovarFGFP Project
## This will be done using the genepi.utils package

#Lib
library("dplyr")
library("genepi.utils")
library("MungeSumstats")
library("gwasvcf")
library("VariantAnnotation")
library("ggplot2")
library("ggpubr")
library("locuscomparer")


#Data path for proj dir
data.path <- "../../../data/PhenoPRS"

#Also load in path to scratch, proj_dir and fgfp data
source("../../parameters/base_dir.R")
source("../../parameters/scratch_path.txt")
fullsspath <- file.path(data.path, "davesRNTsumstats")
sspath_mbg <- file.path(data.path, "mibiogen_ss")

##################
# Load GWAS Data
##################

#Load in the data from the enrichment analysis
load(file.path(data.path, "data_out/Enrichment_SNPs_FGFP_MBG.RData"))

#Want to read this full sumstat data in "G_unclassified_P_Firmicutes_RNT"
fgfp_gwas <- GWAS(dat = file.path(fullsspath, "G_unclassified_P_Firmicutes_RNT_allchr.txt.gz"),
                  map = c(rsid = "rsid", chr = "chromosome", bp = "position", ea = "alleleB", oa = "alleleA", eaf = "all_maf", beta = "frequentist_add_beta_1", se = "frequentist_add_se_1", p = "frequentist_add_pvalue"), 
                  trait = "G_unclassified_P_Firmicutes_RNT", 
                  id    = "fgfpmicro_gwas")

#Lets add in the MiBioGen trait with evidence of genome-wide sig hit with
mbg_gwas <- GWAS(dat = file.path(sspath_mbg, "genus..Ruminococcustorquesgroup.id.14377.summary.txt.gz"),
                 map = c(rsid = "rsID", chr = "chr", bp = "bp", ea = "eff.allele", oa = "ref.allele", beta = "beta", eaf = "eaf", se = "SE", p = "P.weightedSumZ"), 
                 fill = T,
                 trait = "genus..Ruminococcustorquesgroup.id.14377", 
                 id    = "mbg_micro_gwas")


#Now read in the full ggt summary statistic data
if(!file.exists(file.path(data.path, 
                          "phenosumstats/Gamma_glutamyltransferase_levels_UKB_data_field_30730_/ggt_sumstats_vcf_converted.csv.gz"))){

vcffile <- file.path(data.path, 
                     "phenosumstats/Gamma_glutamyltransferase_levels_UKB_data_field_30730_/ebi-a-GCST90014004.vcf.gz")
vcf <- readVcf(vcffile)
vcf_df <- vcf_to_tibble(vcf, id = NULL) 

vcf_df <- vcf_df %>%
  dplyr::select(rsid, seqnames, start, REF, ALT, ES, SE, LP, AF) %>%
  dplyr::rename(
    chr = seqnames,
    pos = start,
    ref = REF,
    alt = ALT,
    beta = ES,
    se = SE,
    lp = LP,
    eaf = AF
  ) %>%
  dplyr::mutate(
    pval = 10^(-lp)
  )

#Output the df as a .gz file
write.csv(vcf_df, 
          file = file.path(data.path, 
                           "phenosumstats/Gamma_glutamyltransferase_levels_UKB_data_field_30730_/ggt_sumstats_vcf_converted.csv.gz"), 
          row.names = FALSE)
}

ggt_gwas <- GWAS(dat = file.path(data.path, 
                                        "phenosumstats/Gamma_glutamyltransferase_levels_UKB_data_field_30730_/ggt_sumstats_vcf_converted.csv.gz"),
                 map = c(
                   rsid = "rsid",
                   chr = "chr",
                   bp = "pos",
                   ea = "alt",
                   oa = "ref",
                   eaf = "eaf",
                   beta = "beta",
                   se = "se",
                   p = "pval"
                 ),
                  trait = "ggt", 
                  id    = "ggt_gwas")

##################
# Manhattan plots
##################

### Now can make stacked Manhattan plots
highlight_snps_fgfp <- fgfp_gwas@rsid[[which.min(fgfp_gwas@p)]]
annotate_snps_fgfp <- fgfp_gwas@rsid[[which.min(fgfp_gwas@p)]]


### Now can make stacked Manhattan plots
highlight_snps_mbg <- mbg_gwas@rsid[[which.min(mbg_gwas@p)]]
annotate_snps_mbg <- mbg_gwas@rsid[[which.min(mbg_gwas@p)]]

#Old code for below for one or two combined manhattan plots


#man_fgfp <- manhattan(fgfp_gwas, 
#                      highlight_snps = highlight_snps,
#                      highlight_win = 250,
#                      annotate_snps = annotate_snps,
#                      hit_table = FALSE,
#                      title = "Manhattan Plot of G_unclassified_P_Firmicutes_RNT from Hughes et al")
#
##Output this plot
#ggsave(
#  filename = file.path(data.path, "manplots_full/Manhattan_FGFP_G_unclassified_P_Firmicutes_RNT.png"),
#  plot = man_fgfp,
#  width = 15,      # optional: set dimensions
#  height = 6,      # optional
#  dpi = 600        # optional: high-resolution output
#)


#Lets also make a man plot for ggt so can combine to miami
#man_ggt <- manhattan(ggt_gwas, 
#                      highlight_snps = highlight_snps,
#                      highlight_win = 250,
#                      annotate_snps = annotate_snps,
#                      hit_table = FALSE,
#                      title = "Manhattan Plot of GGT GWAS")
#
##Output this plot
#ggsave(
#  filename = file.path(data.path, "manplots_full/Manhattan_GGT_GWAS.png"),
#  plot = man_ggt,
#  width = 15,      # optional: set dimensions
#  height = 6,      # optional
#  dpi = 600        # optional: high-resolution output
#)
#
#colours_top    <- c("#67A3D9","#C8E7F5")
#colours_bottom <- c("#F8B7CD","#F6D2E0")
#
#source("Miami2_func.R")
#
#p <- miami2(gwases         = list("top"=fgfp_gwas, "bottom"=ggt_gwas),
#           highlight_snps = list("top"=highlight_snps, "bottom"=highlight_snps),
#           highlight_win  = list("top"=250,"bottom"=250),
#           annotate_snps  = list("top"=annotate_snps, "bottom"=annotate_snps),
#           colours        = list("top"=colours_top, "bottom"=colours_bottom),
#           downsample     = 0.4, 
#           hit_table      = FALSE,
#           title          = "Miami Plot",
#           subtitle       = list("top"="Manhattan Plot of G_unclassified_P_Firmicutes_RNT from Hughes et al", 
#                                 "bottom"= "Manhattan Plot of GGT GWAS"),
#           y_limits = list(top =c(0, 8), bottom=c(0, 100)))
#
##Output this miami plot
#ggsave(
#  filename = file.path(data.path, "manplots_full/Miami_FGFP_GGT.png"),
#  plot = p,
#  width = 10,      # optional: set dimensions
#  height = 6,      # optional
#  dpi = 300        # optional: high-resolution output
#)
#


# Three Manhattan plots combined

# Create upper plot
plot_upper <- manhattan(
  gwas             = mbg_gwas,
  highlight_snps   = highlight_snps_mbg,
  highlight_win    = 250,
  annotate_snps    = annotate_snps_mbg,
  colours          = c("#67A3D9","#C8E7F5"),
  highlight_colour = "#FF7043",
  highlight_shape  = 16,
  sig_line_1       = 5e-8,
  sig_line_2       = 5e-8,
  y_limits         = c(0,10),
  hit_table        = FALSE,
  #max_table_hits   = max_table_hits,
  downsample       = 0.4,
  downsample_pval  = 0.1
) + ggplot2::theme(
  axis.title.x = ggplot2::element_blank(),
  plot.background = ggplot2::element_rect(fill='white', color=NA),
  panel.background = ggplot2::element_rect(fill='white')
) +
  ggplot2::labs(
    subtitle = "GWAS 1: MBG_genus.Ruminococcustorquesgroup.id.14377"   
  )

#Now create middle plot for FGFP
plot_middle <- manhattan(
  gwas             = fgfp_gwas,
  highlight_snps   = highlight_snps_fgfp,
  highlight_win    = 250,
  annotate_snps    = annotate_snps_fgfp,
  colours          = c("#F8B7CD","#F6D2E0"),
  highlight_colour = "#FF7043",
  highlight_shape  = 16,
  sig_line_1       = 5e-8,
  sig_line_2       = 5e-8,
  y_limits         = c(0,10),
  hit_table        = FALSE,
  #max_table_hits   = max_table_hits,
  downsample       = 0.4,
  downsample_pval  = 0.1
) + ggplot2::theme(
  axis.title.x = ggplot2::element_blank(),
  plot.background = ggplot2::element_rect(fill='white', color=NA),
  panel.background = ggplot2::element_rect(fill='white')
) +
  ggplot2::labs(
    subtitle = "GWAS 2: FGFP_G_unclassified_P_Firmicutes"   
  )

# Create lower plot
plot_lower <- manhattan(
  gwas             = ggt_gwas,
  highlight_snps   = c(highlight_snps_mbg, highlight_snps_fgfp),
  highlight_win    = 250,
  annotate_snps    = c(annotate_snps_fgfp, annotate_snps_mbg),
  colours          = c("#8DD3C7", "#C7EAE5"),
  highlight_colour = "#FF7043",
  highlight_shape  = 16,
  sig_line_1       = 5e-8,
  sig_line_2       = 5e-8,
  y_limits         = c(0,130),
  hit_table        = FALSE,
  #max_table_hits   = max_table_hits,
  downsample       = 0.4,
  downsample_pval  = 0.1
)


# Flip bottom plot
plot_lower <- plot_lower +
  ggplot2::scale_y_reverse(limits=rev(c(0,130)), expand=c(0,0)) +
  ggplot2::scale_x_continuous(expand=c(0.01,0.01), position="top") +
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.line.x  = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    axis.text.x  = ggplot2::element_blank(),
    plot.background  = ggplot2::element_rect(fill='white', color=NA),
    panel.background = ggplot2::element_rect(fill='white')
  ) +
  ggplot2::labs(
    subtitle = "GWAS 3: GGT"
  )


full_plot <- plot_upper / plot_middle / plot_lower +
  plot_annotation(
    title = "Miami Plot of Shared Loci across Microbiome Traits and Serum GGT",
    theme = theme(
      plot.title = element_text(face="bold", size=16)
    )
  )

#Output plot
ggsave(
  filename = file.path(data.path, "manplots_full/Miami3_MBG_FGFP_GGT.png"),
  plot = full_plot,
  width = 11,      # optional: set dimensions
  height = 9,      # optional
  dpi = 300        # optional: high-resolution output
)

##################
# Locus plots
##################

#Lets also make locus plots so can get can an idea of beta directional between the two traits of interest
#Lets get the harm_func up and running
source("Harm_func.R")

harm_ggt_GunPFirm <- harmonise_locusRegion(SNP = highlight_snps_fgfp,
                      win = 100,
                      gwas1 = fgfp_gwas,
                      gwas2 = ggt_gwas)

harm_ggt_GRumi <- harmonise_locusRegion(SNP = highlight_snps_mbg,
                                     win = 100,
                                     gwas1 = mbg_gwas,
                                     gwas2 = ggt_gwas)

cleanharm1 <- list(harm_ggt_GunPFirm = harm_ggt_GunPFirm,
                   harm_ggt_GRumi = harm_ggt_GRumi)

#Save this so can work with the data locally
save(cleanharm1, file = file.path(data.path, "manplots_full/Harmonised_Locus1.RData"))
load(file.path(data.path, "manplots_full/Harmonised_Locus1.RData"))


#Now can generate the locus plots
#Make the beta directional plot 
source("LocusBeta_func.R")

plot_data_1 <- data.frame(rsid = cleanharm1$harm_ggt_GunPFirm$SNP, 
                          pval = cleanharm1$harm_ggt_GunPFirm$pval.exposure,
                          beta = cleanharm1$harm_ggt_GunPFirm$beta.exposure
                          )

plot_data_2 <- data.frame(rsid = cleanharm1$harm_ggt_GunPFirm$SNP,
                          pval = cleanharm1$harm_ggt_GunPFirm$pval.outcome,
                          beta = cleanharm1$harm_ggt_GunPFirm$beta.outcome
                          )

#coloc_plots <- locuscomparer::locuscompare(in_fn1 = plot_data_1, 
#                                           in_fn2 = plot_data_2, 
#                                           marker_col1 = "rsid", 
#                                           pval_col1 = "pval", 
#                                           title1 = "G_unclassified_P_Firmicutes_FGFP", 
#                                           marker_col2 = "SNP", 
#                                           pval_col2 = "pval",
#                                           title2 = "rsid", 
#                                           snp = plot_data_1$rsid[which.min(plot_data_1$pval)], 
#                                           population = "EUR", 
#                                           combine = FALSE, 
#                                           legend = TRUE, 
#                                           legend_position = "topright", 
#                                           lz_ylab_linebreak = FALSE)


#Lets also add in the mibiogen trait where there was enrichment for the GGT SNP
beta1 <- betalocuscompare(in_fn1 = plot_data_1, 
                 in_fn2 = plot_data_2, 
                 marker_col1 = "rsid", 
                 pval_col1 = "pval", 
                 beta_col1 = "beta",
                 title1 = "G_unclassified_P_Firmicutes_FGFP", 
                 marker_col2 = "SNP", 
                 pval_col2 = "pval", 
                 beta_col2 = "beta",
                 title2 = "GGT", 
                 snp = plot_data_1$rsid[which.min(plot_data_1$pval)], 
                 population = "EUR", 
                 combine = FALSE, 
                 legend = T, 
                 legend_position = "bottomright")

#Lets output this plot A6 size landscape
ggsave(
  filename = file.path(data.path, "manplots_full/Locus_Beta_GunPFirm_GGT.pdf"),
  plot = beta1,
  width = 5.8,      
  height = 4.1,      
  dpi = 600        
)


#Now do the same for the mbg trait
plot_data_3 <- data.frame(rsid = cleanharm1$harm_ggt_GRumi$SNP, 
                          pval = cleanharm1$harm_ggt_GRumi$pval.exposure,
                          beta = cleanharm1$harm_ggt_GRumi$beta.exposure
)

plot_data_4 <- data.frame(rsid = cleanharm1$harm_ggt_GRumi$SNP,
                          pval = cleanharm1$harm_ggt_GRumi$pval.outcome,
                          beta = cleanharm1$harm_ggt_GRumi$beta.outcome
)


beta2 <- betalocuscompare(in_fn1 = plot_data_3, 
                 in_fn2 = plot_data_4, 
                 marker_col1 = "rsid", 
                 pval_col1 = "pval", 
                 beta_col1 = "beta",
                 title1 = "genus..Ruminococcustorquesgroup_MBG", 
                 marker_col2 = "SNP", 
                 pval_col2 = "pval", 
                 beta_col2 = "beta",
                 title2 = "GGT", 
                 snp = plot_data_3$rsid[which.min(plot_data_3$pval)], 
                 population = "EUR", 
                 combine = FALSE, 
                 legend = T, 
                 legend_position = "topright")

#Now output this plot
ggsave(
  filename = file.path(data.path, "manplots_full/Locus_Beta_GRumi_GGT.pdf"),
  plot = beta2,
  width = 5.8,      
  height = 4.1,      
  dpi = 600        
)

#coloc_plots <- locuscomparer::locuscompare(in_fn1 = plot_data_3, 
#                                           in_fn2 = plot_data_4, 
#                                           marker_col1 = "rsid", 
#                                           pval_col1 = "pval", 
#                                           title1 = "G_unclassified_P_Firmicutes_FGFP", 
#                                           marker_col2 = "SNP", 
#                                           pval_col2 = "pval",
#                                           title2 = "rsid", 
#                                           snp = plot_data_3$rsid[which.min(plot_data_3$pval)], 
#                                           population = "EUR", 
#                                           combine = FALSE, 
#                                           legend = TRUE, 
#                                           legend_position = "topright", 
#                                           lz_ylab_linebreak = FALSE)









#######################################
# Lets do this again for the Lach LOCI
#######################################

#Lets add in the MiBioGen trait with evidence of genome-wide sig hit with
mbg_gwas2 <- GWAS(dat = file.path(sspath_mbg, "family.Lachnospiraceae.id.1987.summary.txt.gz"),
                 map = c(rsid = "rsID", chr = "chr", bp = "bp", ea = "eff.allele", oa = "ref.allele", beta = "beta", eaf = "eaf", se = "SE", p = "P.weightedSumZ"), 
                 fill = T,
                 trait = "family.Lachnospiraceae.id.1987", 
                 id    = "mbg_micro_gwas")

covar_gwas <- list()
for(covar in filter(summary_table_mbg, Micro_Trait == "family.Lachnospiraceae") %>% pull(Covar)){
  #Now read in the full summary statistic data
  sspath <- file.path(data.path, 
                      "phenosumstats", covar)
  
  #File will end in .gz 
  file <- list.files(sspath)[grep(".gz$", list.files(sspath))]
  
  # If length(file) == 1 then move on, otherwise let user choose
  if (length(file) == 1) {
    # Only one file, proceed
    file <- file
  } else {
    # Multiple files detected, prompt user for choice
    cat("Multiple files detected. Please choose one of the following:\n")
    for (i in seq_along(file)) {
      cat(paste0("[", i, "] ", basename(file[i]), "\n"))
    }
    
    selection <- as.integer(readline(prompt = "Enter the number of the file to use: "))
    
    # Check for valid input
    if (is.na(selection) || selection < 1 || selection > length(file)) {
      stop("Invalid selection. Please rerun and choose a valid number.")
    }
    
    file <- file[selection]
  }
  
    
    
    if(covar == "yengo_bmi"){
      map <- c(
        rsid = "SNP",
        chr = "CHR",
        bp = "POS",
        ea = "Tested_Allele",
        oa = "Other_Allele",
        eaf = "Freq_Tested_Allele_in_HRS",
        beta = "BETA",
        se = "SE",
        p = "P"
      )
    }else if(covar == "triglycerides"){
      map <- c(
        rsid = "rsid",
        chr = "chromosome",
        bp = "base_pair_location",
        ea = "effect_allele",
        oa = "other_allele",
        eaf = "effect_allele_frequency",
        beta = "beta",
        se = "standard_error",
        p = "p_value"
      )
    }else{
      map <- c(
        rsid = "hm_rsid",
        chr = "hm_chrom",
        bp = "hm_pos",
        ea = "hm_effect_allele",
        oa = "hm_other_allele",
        eaf = "hm_effect_allele_frequency",
        beta = "hm_beta",
        se = "standard_error",
        p = "p_value"
      )
    }
  
  covar_gwas[[covar]] <- GWAS(dat = file.path(data.path, 
                                   "phenosumstats", covar, file),
                   map = map,
                   trait = covar, 
                   id    = paste0(covar, "_GWAS"),
                   remove = T)
}

#Now get the SNP of inter
SNP4plots <- filter(summary_table_mbg, Micro_Trait == "family.Lachnospiraceae") %>% pull(Micro_rsid) %>% unique()


#Now harmonise all the covariate GWAS for SNP of interest like have done previously
harm_lach_covar <- list()
for(covar in names(covar_gwas)){
  harm_lach_covar[[covar]] <- harmonise_locusRegion(SNP = SNP4plots,
                                                    win = 1000,
                                                    gwas1 = mbg_gwas2,
                                                    gwas2 = covar_gwas[[covar]])
}

#Save this and load locally
save(harm_lach_covar, file = file.path(data.path, "manplots_full/Harmonised_Lach_Covar_Locus.RData"))
load(file.path(data.path, "manplots_full/Harmonised_Lach_Covar_Locus.RData"))

#Now can make locus plots with Pvalue and beta direction for the microbiome trait with the covariates
plots <- list()
for(covar in filter(summary_table_mbg, Micro_Trait == "family.Lachnospiraceae") %>% pull(Covar)){
  plot_data_1 <- data.frame(rsid = harm_lach_covar[[covar]]$SNP, 
                            pval = harm_lach_covar[[covar]]$pval.exposure,
                            beta = harm_lach_covar[[covar]]$beta.exposure
  )
  
  plot_data_2 <- data.frame(rsid = harm_lach_covar[[covar]]$SNP,
                            pval = harm_lach_covar[[covar]]$pval.outcome,
                            beta = harm_lach_covar[[covar]]$beta.outcome
  )
  
  plots[[covar]] <- list()
  
  beta1 <- filter(plot_data_1, rsid == plot_data_1$rsid[which.min(plot_data_1$pval)]) %>% pull(beta)
  beta2 <- filter(plot_data_2, rsid == plot_data_1$rsid[which.min(plot_data_1$pval)]) %>% pull(beta)
  
  #Multiple together to get sign
  sign_beta <- beta1 * beta2
  
  if(sign_beta > 0){
    legpos = "bottomright"
  } else{
    legpos = "topright"
  }
  
  #Get the clean name from the covar
  covar_clean <- filter(summary_table_mbg, Covar == covar) %>% pull(covar_clean) %>% unique()
  
  plots[[covar]][["beta"]] <- betalocuscompare(in_fn1 = plot_data_1, 
                            in_fn2 = plot_data_2, 
                            marker_col1 = "rsid", 
                            pval_col1 = "pval", 
                            beta_col1 = "beta",
                            title1 = "family.Lachnospiraceae_MBG", 
                            marker_col2 = "SNP", 
                            pval_col2 = "pval", 
                            beta_col2 = "beta",
                            title2 = covar_clean, 
                            snp = plot_data_1$rsid[which.min(plot_data_1$pval)], 
                            population = "EUR", 
                            combine = FALSE, 
                            legend = T, 
                            legend_position = legpos)
  
  #Now do the locus plot with p-values
  plots[[covar]][["pval"]] <- locuscomparer::locuscompare(in_fn1 = plot_data_1, 
                                             in_fn2 = plot_data_2, 
                                             marker_col1 = "rsid", 
                                             pval_col1 = "pval", 
                                             title1 = "family.Lachnospiraceae_MBG", 
                                             marker_col2 = "SNP", 
                                             pval_col2 = "pval",
                                             title2 = covar_clean, 
                                             snp = plot_data_1$rsid[which.min(plot_data_1$pval)], 
                                             population = "EUR", 
                                             combine = FALSE, 
                                             legend = FALSE, 
                                             legend_position = "topright", 
                                             lz_ylab_linebreak = FALSE)

  
}

#Patchwork beta plots into a grid
beta_grid <- (plots[[1]][["beta"]] + plots[[2]][["beta"]]) /
  (plots[[3]][["beta"]] + plots[[4]][["beta"]])

p_grid <- (plots[[1]][["pval"]][[2]] /
  plots[[1]][["pval"]][[3]] /
  plots[[2]][["pval"]][[3]] /
  plots[[3]][["pval"]][[3]] /
  plots[[4]][["pval"]][[3]])

