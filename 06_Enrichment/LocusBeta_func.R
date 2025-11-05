#' Read association summary statistics from file and append column.
#' The file must contain 2 columns: markers, i.e SNPs, and p-value.
#' The marker column should contain SNP rsIDs.
#'
#' @param in_fn (string) Path to the input file.
#' @param marker_col (string, optional) Name of the marker column. Default: 'rsid'.
#' @param pval_col (string, optional) Name of the p-value column. Default: 'pval'.
#' @examples
#' in_fn = system.file('extdata', 'gwas.tsv', package = 'locuscomparer')
#' d1 = read_metal(in_fn, marker_col = 'rsid', pval_col = 'pval')
#' @export
#' 
#' 
#' 
#' 

#Essential packages
library(ggrepel)
library(cowplot)

read_metal=function(in_fn,marker_col='rsid',pval_col='pval',beta_col='beta'){
  # message('Reading ', in_fn)
  
  if (is.character(in_fn)){
    
    d = read.table(in_fn, header = TRUE, stringsAsFactors = FALSE)
    colnames(d)[which(colnames(d) == marker_col)] = 'rsid'
    colnames(d)[which(colnames(d) == pval_col)] = 'pval'
    
  } else if (is.data.frame(in_fn)){
    
    d = in_fn
    
  } else {
    
    stop('The argument "in_fn" must be a string or a data.frame')
    
  }
  
  d$logp = -log10(d$pval)
  return(d[,c('rsid','pval','logp','beta')])
}

#' Append two columns, chromosome (chr) and position (pos), to the input data.frame.
#'
#' @param x (data.frame) Input data.frame.
#' @param genome (string, optional) Genome assembly, either 'hg19' or 'hg38'. Default: 'hg19'.
#' @examples
#' in_fn = system.file('extdata', 'gwas.tsv', package = 'locuscomparer')
#' d1 = read_metal(in_fn, marker_col = 'rsid', pval_col = 'pval')
#' get_position(d1, genome)
#' @export
get_position=function(x, genome = c('hg19','hg38')){
  
  data(config)
  on.exit(rm(config))
  
  conn = RMySQL::dbConnect(RMySQL::MySQL(),"locuscompare",config$b,config$c,config$a)
  on.exit(RMySQL::dbDisconnect(conn))
  
  stopifnot('rsid' %in% colnames(x))
  
  genome = match.arg(genome)
  
  cmd = sprintf("select rsid, chr, pos from tkg_p3v5a_%s where rsid in ('%s')",genome,paste0(x$rsid,collapse="','"))
  res = DBI::dbGetQuery(conn = conn, statement = cmd)
  y=merge(x,res,by='rsid')
  return(y)
}


#' Retrive SNP pairwise LD from database.
#' SNP pairwise lD are calculated based on 1000 Genomes Project Phase 3 version 5.
#' For storage-efficiency, the output will only include SNPs with r2 > 0.2 with the
#' input SNP.
#' @param chr (string) Chromosome name. e.g. '22'. Notice that the name should not contain 'chr'.
#' @param snp (string) SNP rsID.
#' @param population (string) One of the 5 popuations from 1000 Genomes: 'AFR', 'AMR', 'EAS', 'EUR', and 'SAS'.
#' @examples
#' retrieve_LD('6', 'rs9349379', 'AFR')
#'
#' @export
retrieve_LD = function(chr,snp,population){
  data(config)
  on.exit(rm(config))
  
  conn = RMySQL::dbConnect(RMySQL::MySQL(),"locuscompare",config$b,config$c,config$a)
  on.exit(RMySQL::dbDisconnect(conn))
  
  res1 = DBI::dbGetQuery(
    conn = conn,
    statement = sprintf(
      "select SNP_A, SNP_B, R2
            from tkg_p3v5a_ld_chr%s_%s
            where SNP_A = '%s';",
      chr,
      population,
      snp
    )
  )
  
  res2 = DBI::dbGetQuery(
    conn = conn,
    statement = sprintf(
      "select SNP_B as SNP_A, SNP_A as SNP_B, R2
            from tkg_p3v5a_ld_chr%s_%s
            where SNP_B = '%s';",
      chr,
      population,
      snp
    )
  )
  
  res = rbind(res1,res2)
  return(res)
}

#' Get the lead SNP from the list of SNPs in input data.frame
#' The lead SNP is defined as the SNP with the lowest sum of p-values
#' from the two studies.
#' @param merged (data.frame) Input data.frame, which is a result by merging two
#' association studies.
#' @param snp (string, optional) SNP rsID. If NULL, the function will select the
#' lead SNP based on the sum of p-values from the two studies. If an rsID is supplied,
#' the function will simply return the rsID.
#' @examples
#' # Select the lead SNP
#' in_fn_1 = system.file('extdata', 'gwas.tsv', package = 'locuscomparer')
#' d1 = read_metal(in_fn_1, marker_col = 'rsid', pval_col = 'pval')
#' in_fn_2 = system.file('extdata', 'gwas.tsv', package = 'locuscomparer')
#' d1 = read_metal(in_fn_2, marker_col = 'rsid', pval_col = 'pval')
#' merged = merge(d1, d2, by = "rsid", suffixes = c("1", "2"), all = FALSE)
#' get_lead_snp(merged)
#' @export
get_lead_snp = function(merged, snp = NULL){
  if (is.null(snp)) {
    snp = merged[which.min(merged$pval1 + merged$pval2), 'rsid']
  }
  else {
    if (!snp %in% merged$rsid) {
      stop(sprintf("%s not found in the intersection of in_fn1 and in_fn2.", snp))
    }
  }
  return(as.character(snp))
}

#' Assign color to each SNP according to LD.
#' @param rsid (character vector) A vector of rsIDs on which to assign color.
#' @param snp (string) rsID for lead SNP. This SNP will be colored purple.
#' Other SNPs will be assigned color based on their LD with the lead SNP.
#' @param ld (data.frame) The output from `retrieve_LD()`.
#' @examples
#' # the data.frame merged comes from the example for `get_lead_snp()`.
#' # the data.frame ld comes from the example for `retrieve_LD()`.
#' color = assign_color(rsid = merged$rsid, snp = 'rs9349379', ld)
#' @export
assign_color=function(rsid,snp,ld){
  
  ld = ld[ld$SNP_A==snp,]
  ld$color = as.character(cut(ld$R2,breaks=c(0,0.2,0.4,0.6,0.8,1), labels=c('blue4','skyblue','darkgreen','orange','red'), include.lowest=TRUE))
  
  color = data.frame(rsid, stringsAsFactors = FALSE)
  color = merge(color, ld[, c('SNP_B', 'color')], by.x = 'rsid', by.y = 'SNP_B', all.x = TRUE)
  color[is.na(color$color),'color'] = 'blue4'
  if (snp %in% color$rsid){
    color[rsid == snp,'color'] = 'purple'
  } else {
    color = rbind(color, data.frame(rsid = snp, color = 'purple'))
  }
  
  res = color$color
  names(res) = color$rsid
  
  return(res)
}


#' Add a column of SNP labels to input data.frame
#' @param merged (data.frame) Input data.frame, which is a result by merging two
#' association studies. See the example under `get_lead_snp()` for generation of
#' such data.frame.
#' @param snp (character vector) A vector of SNP rsIDs. If only labeling one SNP,
#' this can also be a single string.
#' @examples
#' # The data.frame merged comes from the example for `get_lead_snp()`.
#' merged = add_label(merged, 'rs9349379')
add_label = function(merged, snp){
  merged$label = ifelse(merged$rsid %in% snp, merged$rsid, '')
  return(merged)
}


#' Make a scatter plot (called the LocusCompare plot).
#' Each axis of the LocusCompare plot represent the -log10(p-value) from
#' an association study. Each point thus represent a SNP. By default, the lead SNP
#' is a purple diamond, whereas the other SNPs are colored according to
#' their LD with the lead SNP.
#' @import ggplot2
#' @import cowplot
#' @param merged (data.frame) An input data.frame which has the following
#' columns: rsid, pval1 (p-value for study 1), logp1 (p-value for study 2),
#' logp1 (log p-value for study 1), logp2 (log p-value for study 2), chr, pos.
#' See the example for `get_lead_snp()` on how to generate this data.frame.
#' @param title1 (string) The title for the x-axis.
#' @param title2 (string) The title for the y-axis.
#' @param color (data.frame) The output from `assign_color()`.
#' @param shape (data.frame) Specification of the shape of each SNP. See example blow on how to generate this data.frame.
#' @param size (data.frame) Specification of the size of each SNP. See example below on how to generate this data.frame.
#' @param legend (boolean) Whether to include the legend.
#' @param legend_position (string, optional) Either 'bottomright','topright', or 'topleft'. Default: 'bottomright'.
#' @examples
#' # The data.frame `merged` comes from the example of `add_label()`.
#' # The data.frame `color` comes from the example of `assign_color()`.
#' snp = 'rs9349379'
#' shape = ifelse(merged$rsid == snp, 23, 21)
#' names(shape) = merged$rsid
#' size = ifelse(merged$rsid == snp, 3, 2)
#' names(size) = merged$rsid
#' make_scatterplot(merged, title1 = 'GWAS', title2 = 'eQTL', color, shape, size)
#' @export

#Here will be a scatterplot of the betas rather than the P-values 
make_scatterplot_beta <- function(merged, title1, title2, color, shape, size, 
                                  legend = TRUE, legend_position = c('bottomright','topright','topleft')) {
  
  # Ensure lead SNP label exists
  if(!"label" %in% colnames(merged)) merged$label <- ""
  
  # Identify SNPs in LD (not "blue4") OR the lead SNP itself
  high_ld_snps <- merged$rsid[color[merged$rsid] != "blue4" | merged$label != ""]
  
  # Filter merged to only those SNPs
  merged <- merged[merged$rsid %in% high_ld_snps, ]
  
  # Recalculate model for high-LD SNPs
  if(length(high_ld_snps) > 1){
    lm_fit <- lm(beta2 ~ beta1, data = merged)
    slope <- coef(lm_fit)[2]
    intercept <- coef(lm_fit)[1]
  } else {
    slope <- 0
    intercept <- 0
  }
  
  # Build plot
  p <- ggplot(merged, aes(beta1, beta2)) +
    geom_point(aes(fill = rsid, size = rsid, shape = rsid), alpha = 0.8) +
    geom_point(data = merged[merged$label != "", ],
               aes(beta1, beta2, fill = rsid, size = rsid, shape = rsid)) +
    geom_abline(intercept = intercept, slope = slope, color = "grey50", linetype = "dashed") +
    ggrepel::geom_label_repel(
      data = merged[merged$label != "", ],
      aes(label = label),
      size = 4,
      fontface = "bold",
      fill = "white",
      color = "black",
      box.padding = 0.25,
      label.size = 0.35,
      segment.color = "grey50"
    ) +
    xlab(bquote(.(title1) ~ 'Beta')) +
    ylab(bquote(.(title2) ~ 'Beta')) +
    scale_fill_manual(values = color, guide = "none") +
    scale_shape_manual(values = shape, guide = "none") +
    scale_size_manual(values = size, guide = "none") +
    theme_classic()
  
  # Optional legend
  if (legend) {
    legend_position <- match.arg(legend_position)
    if (legend_position == 'bottomright'){
      legend_box <- data.frame(x = 0.8, y = seq(0.4, 0.2, -0.05))
    } else if (legend_position == 'topright'){
      legend_box <- data.frame(x = 0.8, y = seq(0.8, 0.6, -0.05))
    } else {
      legend_box <- data.frame(x = 0.2, y = seq(0.8, 0.6, -0.05))
    }
    
    p <- ggdraw(p) +
      geom_rect(data = legend_box,
                aes(xmin = x, xmax = x + 0.05, ymin = y, ymax = y + 0.05),
                color = "black",
                fill = rev(c("blue4", "skyblue", "darkgreen", "orange", "red"))) +
      draw_label("0.8", x = legend_box$x[1] + 0.05, y = legend_box$y[1], hjust = -0.3, size = 10) +
      draw_label("0.6", x = legend_box$x[2] + 0.05, y = legend_box$y[2], hjust = -0.3, size = 10) +
      draw_label("0.4", x = legend_box$x[3] + 0.05, y = legend_box$y[3], hjust = -0.3, size = 10) +
      draw_label("0.2", x = legend_box$x[4] + 0.05, y = legend_box$y[4], hjust = -0.3, size = 10) +
      draw_label(parse(text = "r^2"), x = legend_box$x[1] + 0.05, y = legend_box$y[1], vjust = -2, size = 10)
  }
  
  return(p)
}



betalocuscompare <- function(in_fn1, in_fn2, marker_col1 = "rsid", beta_col1 = "beta", pval_col1 = "pval",
                             title1 = "eQTL", marker_col2 = "rsid", pval_col2 = "pval", beta_col2 = "beta", title2 = "GWAS",
                             snp = NULL, population = "EUR", combine = TRUE, legend = TRUE,
                             legend_position = c('bottomright','topright','topleft'),
                             genome = c('hg19','hg38')) {
  
  d1 <- read_metal(in_fn1, marker_col1, pval_col1, beta_col1)
  d2 <- read_metal(in_fn2, marker_col2, pval_col2, beta_col2)  # fixed typo
  
  merged <- merge(d1, d2, by = "rsid", suffixes = c("1", "2"), all = FALSE)
  genome <- match.arg(genome)
  merged <- get_position(merged, genome)
  
  chr <- unique(merged$chr)
  if (length(chr) != 1) stop('There must be one and only one chromosome.')
  
  snp <- get_lead_snp(merged, snp)
  ld <- retrieve_LD(chr, snp, population)
  
  # Prepare plot aesthetics
  color <- assign_color(merged$rsid, snp, ld)
  merged <- add_label(merged, snp)
  shape <- ifelse(merged$rsid == snp, 23, 21)
  names(shape) <- merged$rsid
  size <- ifelse(merged$rsid == snp, 3, 2)
  names(size) <- merged$rsid
  
  # Generate scatterplot of betas
  p <- make_scatterplot_beta(merged, title1, title2, color, shape, size, legend, legend_position)
  
  return(p)
}
