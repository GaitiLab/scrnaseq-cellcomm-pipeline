options(repos = "https://CRAN.R-project.org")

# 1). Install 'pak' (fast package installer)
message("Installing 'pak'...")
install.packages("pak")

# 2). Install 'scrnaseq-cellcomm'
message("Installing 'scrnaseq-cellcomm'...")
pak::pkg_install("GaitiLab/scrnaseq-cellcomm")

message("Finished!")
