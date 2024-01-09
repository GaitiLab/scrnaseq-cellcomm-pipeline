# 1) create conda environment cpdb
mamba env create -f env/cpdb.yml

# 2) create conda environment cell2cell
mamba env create -f env/cell2cell.yml

# 3) create conda environment main scripts
mamba env create -f env/cci.yml

mamba activate cci
Rscript "env/install.R"
