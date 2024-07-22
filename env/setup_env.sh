# 1) create conda environment cpdb
mamba env create -f cpdb.yml

# 2) create conda environment cell2cell
mamba env create -f cell2cell.yml

# 3) create conda environment main scripts
mamba env create -f cci.yml && source "$HOME/miniforge3/bin/activate" "cci" && Rscript "install.R"

echo "Finished!"