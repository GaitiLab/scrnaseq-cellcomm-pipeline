#!/usr/bin/env bash

# ---- Copy database files ---- #

# Files to be copied
source_interactions_db_dir="001_data_local/interactions_db_v2"
cellchat_db="${source_interactions_db_dir}/cellchat_db.rds"
liana_db="${source_interactions_db_dir}/liana_db.rds"
liana_db_csv="${source_interactions_db_dir}/cell2cell.csv"
ref_db="${source_interactions_db_dir}/ref_db.rds"
cellphone_db=$(find "${source_interactions_db_dir}/cellphonedb_custom" -name "cellphonedb_*.zip")

# Create target directory if not existing
target_interactions_db_dir="data/interactions_db"
mkdir -p ${target_interactions_db_dir}

# Copy files
gcp -f "${cellchat_db}" "${target_interactions_db_dir}/cellchat_db.rds"
gcp -f "${liana_db}" "${target_interactions_db_dir}/liana_db.rds"
gcp -f "${liana_db_csv}" "${target_interactions_db_dir}/cell2cell.csv"
gcp -f  "${ref_db}" "${target_interactions_db_dir}/ref_db.rds"
gcp -f "${cellphone_db}" "${target_interactions_db_dir}/cellphondb.zip"

if [ -e .commit ]
    then
    rm .commit
    git add "${target_interactions_db_dir}"
    git commit --amend -C HEAD --no-verify
fi
exit