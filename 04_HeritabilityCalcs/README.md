## 04_HeritabilityCalcs Contents

This folder contains scripts for analyses involving heritability estimation

* 04.1 Contains scripts to edit the bgen genetic files to filter on only those SNPs we are interested in for the analyses and convert the file format so can be analysed by GREML. The GRM is also created here

* 04.2 This folder contains scripts to extract our heritable MTs, residualise them against covariate PRSs and calculate heritability estimates for the residualised traits

### TO NOTE

* When converting to hard-call genotypes (what happens when you convert to plink format) it seems to be best practice to filter on INFO score and genotype prob score (e.g. >0.3 and >0.9 respectively). In Dave's FGFP methodology for h2 estimation it states filtering has only been done on INFO score > 0.9, I have done the latter as we would like consistency in h2 estimation, but it may be better practice to do the former.

* The samples in bgen files who were flagged for relatedness in Daves GRM have not been removed, instead I have calculated my own GRM at 0.25 cut off as stated in Dave's methodology. There seems to be differences in the individuals removed from Dave's estimation and mine which cause minor differences in the h2 estimation. This has not yet been formally investigated. 
