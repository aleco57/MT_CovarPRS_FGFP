#Do with patchwork instead
library("patchwork")

#Updated Miami plot
hit_table <- function(gwas, n) {
  # base columns and BETA and SE if provided
  cols <- c(c("rsid","p"), names(gwas)[names(gwas) %in% c("beta","se")])
  table_data <- utils::head(gwas[order(p), cols, with=FALSE], n=n)
  table_data <- table_data[, lapply(.SD, signif, digits=2), by=rsid]
  table_data[, rsid := strtrim(rsid, 23)]
  table <- gridExtra::tableGrob(table_data, rows=NULL)
  return(table)
}


# Updated Miami plot function with consistent white backgrounds
miami2 <- function(gwases,
                   highlight_snps   = list("top"=NULL, "bottom"=NULL),
                   highlight_win    = list("top"=100, "bottom"=100),
                   annotate_snps    = list("top"=NULL,"bottom"=NULL),
                   colours          = list("top"=c("#d9d9d9","#bfbfbf"),"bottom"=c("#bfbfbf","#d9d9d9")),
                   highlight_colour = list("top"="#e15758","bottom"="#4f79a7"),
                   highlight_shape  = list("top"=16,"bottom"=16),
                   sig_line_1       = list("top"=5e-8,"bottom"=5e-8),
                   sig_line_2       = list("top"=NULL,"bottom"=NULL),
                   y_limits         = list("top"=c(NULL,NULL),"bottom"=c(NULL,NULL)),
                   title            = NULL,
                   subtitle         = list("top"=NULL,"bottom"=NULL),
                   base_text_size   = 14,
                   hit_table        = FALSE,
                   max_table_hits   = 10,
                   downsample       = 0.1,
                   downsample_pval  = 0.1) {
  
  # Checks
  stopifnot("Ensure inputs are lists of length two" = all(sapply(
    list(gwases, highlight_snps, highlight_win, annotate_snps, colours, highlight_colour, 
         highlight_shape, sig_line_1, sig_line_2, y_limits),
    function(x) is.list(x) & length(x) == 2
  )))
  
  # Create upper plot
  plot_upper <- manhattan(
    gwas             = gwases[[1]],
    highlight_snps   = highlight_snps[[1]],
    highlight_win    = highlight_win[[1]],
    annotate_snps    = annotate_snps[[1]],
    colours          = colours[[1]],
    highlight_colour = highlight_colour[[1]],
    highlight_shape  = highlight_shape[[1]],
    sig_line_1       = sig_line_1[[1]],
    sig_line_2       = sig_line_2[[1]],
    y_limits         = y_limits[[1]],
    hit_table        = FALSE,
    max_table_hits   = max_table_hits,
    downsample       = downsample,
    downsample_pval  = downsample_pval
  ) + ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    plot.background = ggplot2::element_rect(fill='white', color=NA),
    panel.background = ggplot2::element_rect(fill='white')
  ) +
   ggplot2::labs(
     subtitle = subtitle[[1]]   
    )
  
  # Create lower plot
  plot_lower <- manhattan(
    gwas             = gwases[[2]],
    highlight_snps   = highlight_snps[[2]],
    highlight_win    = highlight_win[[2]],
    annotate_snps    = annotate_snps[[2]],
    colours          = colours[[2]],
    highlight_colour = highlight_colour[[2]],
    highlight_shape  = highlight_shape[[2]],
    sig_line_1       = sig_line_1[[2]],
    sig_line_2       = sig_line_2[[2]],
    y_limits         = y_limits[[2]],
    hit_table        = FALSE,
    max_table_hits   = max_table_hits,
    downsample       = downsample,
    downsample_pval  = downsample_pval
  )
  
  # Recalculate flipped y-limits for bottom plot
  if(is.null(y_limits[[2]])) {
    max_p <- max(-log10(gwases[[2]]@p), na.rm=TRUE)
    if(!is.null(sig_line_1[[2]])) max_p <- max(max_p, -log10(sig_line_1[[2]]), na.rm=TRUE)
    if(!is.null(sig_line_2[[2]])) max_p <- max(max_p, -log10(sig_line_2[[2]]), na.rm=TRUE)
    y_limits[[2]] <- c(ceiling(max_p), 0)
  } else {
    y_limits[[2]] <- rev(y_limits[[2]])
  }
  
  # Flip bottom plot
  plot_lower <- plot_lower +
    ggplot2::scale_y_reverse(limits=y_limits[[2]], expand=c(0,0)) +
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
      subtitle = subtitle[[2]]   
    )
  
  # Combine plots
  # Combine top and bottom plots (without tables)
  #plot <- ggpubr::ggarrange(
  #  plot_upper, 
  #  plot_lower, 
  #  ncol = 1, 
  #  heights = c(1,1), 
  #  bg = "white"  # ensures white background
  #)
  
  plot <- plot_upper / plot_lower +
    plot_annotation(
      title = "Miami Plot",
      theme = theme(
        plot.title = element_text(face="bold", size=16)
      )
    )
  
  # Add top subtitle, main title, and bottom subtitle using annotate_figure
  #plot <- ggpubr::annotate_figure(
  #  plot,
  #  top = ggpubr::ggarrange(
  #    ggpubr::text_grob(subtitle[[1]], size = base_text_size),
  #    ggpubr::text_grob(title, face = "bold", size = base_text_size + 2),
  #    ncol = 1, heights = c(1,2)
  #  ),
  #  bottom = ggpubr::text_grob(subtitle[[2]], size = base_text_size)
  #)
  
  
  return(plot)
}



### Also make a miami3 funciton which can take three GWAS objects of interest to generate the plot
miami3 <- function(gwases,
                   highlight_snps   = list("top1"=NULL, "top2"=NULL, "bottom"=NULL),
                   highlight_win    = list("top1"=100, "top2"=100, "bottom"=100),
                   annotate_snps    = list("top1"=NULL,"top2"=NULL,"bottom"=NULL),
                   colours          = list("top1"=c("#d9d9d9","#bfbfbf"),
                                           "top2"=c("#cccccc","#999999"),
                                           "bottom"=c("#bfbfbf","#d9d9d9")),
                   highlight_colour = list("top1"="#e15758",  # red
                                           "top2"="#57b28f",  # teal
                                           "bottom"=c("#e15758","#57b28f")), # both colours
                   highlight_shape  = list("top1"=16,"top2"=16,"bottom"=16),
                   sig_line_1       = list("top1"=5e-8,"top2"=5e-8,"bottom"=5e-8),
                   y_limits         = list("top1"=c(NULL,NULL),"top2"=c(NULL,NULL),"bottom"=c(NULL,NULL)),
                   title            = "Triple Miami Plot",
                   subtitle         = list("top1"=NULL,"top2"=NULL,"bottom"=NULL),
                   base_text_size   = 14,
                   downsample       = 0.1,
                   downsample_pval  = 0.1) {
  
  stopifnot(length(gwases) == 3)
  
  # ---- Helper to make a single Manhattan ----
  make_manh <- function(gwas, which, flipped = FALSE) {
    p <- manhattan(
      gwas             = gwas,
      highlight_snps   = highlight_snps[[which]],
      highlight_win    = highlight_win[[which]],
      annotate_snps    = annotate_snps[[which]],
      colours          = colours[[which]],
      highlight_colour = highlight_colour[[which]],
      highlight_shape  = highlight_shape[[which]],
      sig_line_1       = sig_line_1[[which]],
      y_limits         = y_limits[[which]],
      hit_table        = FALSE,
      downsample       = downsample,
      downsample_pval  = downsample_pval
    ) + 
      theme(
        plot.background = element_rect(fill='white', color=NA),
        panel.background = element_rect(fill='white'),
        axis.title.x = element_blank()
      ) +
      labs(subtitle = subtitle[[which]])
    
    if (flipped) {
      # reverse y-axis for bottom plot
      max_p <- max(-log10(gwas$p), na.rm=TRUE)
      ylim <- if(is.null(y_limits[[which]])) c(ceiling(max_p), 0) else rev(y_limits[[which]])
      p <- p +
        scale_y_reverse(limits = ylim, expand = c(0,0)) +
        scale_x_continuous(expand = c(0.01,0.01), position = "top") +
        theme(
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank()
        )
    }
    return(p)
  }
  
  # ---- Create plots ----
  plot_top1 <- make_manh(gwases[[1]], "top1")
  plot_top2 <- make_manh(gwases[[2]], "top2")
  plot_bottom <- make_manh(gwases[[3]], "bottom", flipped = TRUE)
  
  # ---- Combine ----
  plot <- (plot_top1 / plot_top2 / plot_bottom) +
    plot_annotation(
      title = title,
      theme = theme(plot.title = element_text(face="bold", size=16))
    )
  
  return(plot)
}


