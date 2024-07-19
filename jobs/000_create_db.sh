#!/usr/bin/env bash
source "$HOME/miniforge3/bin/activate" "cpdb"

python "${1}/Python/017_update_cellphonedb.py" \
    --input_dir $2


