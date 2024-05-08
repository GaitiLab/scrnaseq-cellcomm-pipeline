# Visualization
devtools::install_github("jmw86069/colorjam", upgrade = "never")

# CellChat (order is important, presto necessary because of bug w/ CellChat installation)
devtools::install_github("immunogenomics/presto", upgrade = "never")
devtools::install_github("jinworks/CellChat", upgrade = "never")

# LIANA
devtools::install_github("saezlab/liana", upgrade = "never")


devtools::install_github("GaitiLab/GaitiLabUtils")
devtools::install_github("GaitiLab/GBMutils")

# Install scalop for signature
# Install first according to https://github.com/jlaffy/scalop/issues/5#issuecomment-1640314383
BiocManager::install("Homo.sapiens")
# devtools::install_github("jlaffy/scalop") # unable to install
devtools::install_github("jpmam1/scalop", upgrade = "never")
# use-fork as suggested here: https://stackoverflow.com/questions/76387189/fatal-error-relating-to-include-s-h-when-installing-r-scalop-package/76387609#76387609
