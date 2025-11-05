#### Function for genepi.utils which can select a region of interest from sumstats for locus plots

#What is the SNP of interest?
harmonise_locusRegion <- function(SNP = SNP,
                       win = 100,
                       gwas1 = gwas1,
                       gwas2 = gwas2){
#Testing
if(FALSE){
  SNP = "rs2304779"
  win = 100
  gwas1 = fgfp_gwas
  gwas2 = ggt_gwas
  title1 = "G_unclassified_P_Firmicutes_FGFP"
  title2 = "GGT"
}

  #Convert to df so can work with easier
  gwas1 <- genepi.utils::as.data.table(gwas1)
  gwas2 <- genepi.utils::as.data.table(gwas2)
  
  bp_  <- gwas1$bp[gwas1$rsid == SNP]
  chr_ <- gwas1$chr[gwas1$rsid == SNP]
  
  #Check this is the same for gwas2, print an error if not
  bp2_  <- gwas2$bp[gwas2$rsid == SNP]
  chr2_ <- gwas2$chr[gwas2$rsid == SNP]
  
  if(length(bp_) < 1 | length(chr_) < 0){
    #stop and check if supplied the right SNP
    stop("Error: Macth SNP not availble in GWAS1, please check SNP ID")
  }else if(length(bp2_) < 1 | length(chr2_) < 0){
    print("Warning: Macth SNP not availble in GWAS2, will use GWAS1 ONLY")
    rsidstomatch <- gwas1[bp > (bp_ - win*1000) & bp < (bp_ + win*1000) & chr == chr_, "rsid"]
    subset1 <- gwas1[rsid %in% rsidstomatch$rsid, ]
    subset2 <- gwas2[rsid %in% rsidstomatch$rsid, ]
  }else if(bp_ != bp2_ | chr_ != chr2_){
    print("SNP position does not match between gwas1 and gwas2, are these the correct genome build?\n")
    print("Will match of RSIDs instead\n")
    
    #Match off rsids from subset 1
    rsidstomatch <- gwas1[bp > (bp_ - win*1000) & bp < (bp_ + win*1000) & chr == chr_, "rsid"]
    subset1 <- gwas1[rsid %in% rsidstomatch$rsid, ]
    subset2 <- gwas2[rsid %in% rsidstomatch$rsid, ]
  }else{
  
  subset1 <- gwas1[bp > (bp_ - win*1000) & bp < (bp_ + win*1000) & chr == chr_,]
  subset2 <- gwas2[bp > (bp_ - win*1000) & bp < (bp_ + win*1000) & chr == chr_,]
  
  }
  #Now lets harmonise these using the TwoSampleMR harmonise_data function
  harm <- TwoSampleMR::harmonise_data(
    TwoSampleMR::format_data(data.frame(subset1), type = "exposure", phenotype_col = "id", snp_col = "rsid", beta_col = "beta", se_col = "se",
                eaf_col = "eaf", effect_allele_col = "ea", other_allele_col = "oa", pval_col = "p",
                chr_col = "chr", pos_col = "bp"),
    TwoSampleMR::format_data(data.frame(subset2), type = "outcome", phenotype_col = "id", snp_col = "rsid", beta_col = "beta", se_col = "se",
                             eaf_col = "eaf", effect_allele_col = "ea", other_allele_col = "oa", pval_col = "p",
                             chr_col = "chr", pos_col = "bp")
  )
    
    
  #Now can return the harmonised locus region
  return(harm)
}



