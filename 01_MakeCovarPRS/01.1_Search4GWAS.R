library(ieugwasr)
library(dplyr)

#Collect all metadata in IEUGWAS data
info <- gwasinfo()

#Search for GWAS iteratively using grep and filtering to only european and biggest sample size
search <- "Ulcerative colitis"
search_df <- info[grep(tolower(search), tolower(info$trait)), ] %>% as.data.frame()

search_df <- search_df[search_df$population == "European", ]

# Order the results by the sample_size column in decreasing order
search_df <- search_df[order(search_df$sample_size, decreasing = TRUE), ]

View(search_df)
