# Rscript to extract all SNPs 
#This needs to be run in bp directory

#Note the paths have been removed, to see the paths look on bp script folder

library(dplyr)
library(data.table)

file_dir <- ""

all_files <- list.files(file_dir) %>%
  #Remove BetaDiv
  grep("BetaDiv", ., invert = TRUE, value = TRUE)

#Lets starts with RNT
rnt_files <- all_files[grep("RNT", all_files)]

out <- list()

for(file in rnt_files) {
  file_path <- file.path(file_dir, file)
  
  #Name
  name <- gsub("_allchr.txt.gz", "", file)
  
  # Read the file
  out[[name]] <- data.table::fread(file = file.path(file_dir, file)) %>% 
    filter(frequentist_add_pvalue < 1e-4)
  
  cat("Processed file:", file, "\n")
}

#Save to the scratch
save(out, 
     file = file.path(""))
  
  
#Now extract HB (probably won't use but good to do now)
hb_files <- all_files[grep("HB", all_files)]
for(file in hb_files) {
  file_path <- file.path(file_dir, file)
  
  #Name
  name <- gsub("_allchr.txt.gz", "", file)
  
  # Read the file
  out[[name]] <- data.table::fread(file = file.path(file_dir, file)) %>% 
    filter(frequentist_add_pvalue < 1e-4)
  
  cat("Processed file:", file, "\n")
}

#Save out again
save(out, 
     file = file.path(""))