############
# TO NOTE !
############

#1.
####
#When converting to hard-call genotypes (what happens when you convert to plink format) it seems to be best practice to filter on INFO score and genotype prob score (e.g. >0.3 and >0.9 respectively)
#In Dave's methodology it states filtering has only been done on INFO score > 0.9, I have done the same here as we would like consistency in h2 estimation, but it may be better practice to do the former

#2.
####
# The samples in bgen files who were flagged for relatedness have not been removed, this may cause potential issues in the analysis - check with Dave if this should be removed
# Note this might not be an issue for calculating GRM as this is exactly what it is trying to calculate
