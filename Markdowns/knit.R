# Code to knit .Rmd in a different directory
rmarkdown::render("h2_mts_res.Rmd", output_dir = "../../../data/PhenoPRS/markdowns/h2MTs_rmd")

rmarkdown::render("Results_update.Rmd", output_dir = "../../../data/PhenoPRS/markdowns/h2MTs_rmd")


#For all the plots
rmarkdown::render("CovarPRS_plots.Rmd", output_dir = "../../../data/PhenoPRS/markdowns/plots_rmd")
