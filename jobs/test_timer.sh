#!/usr/bin/env bash
work_dir="/Users/joankant/Desktop/gaitigroup/Users/Joan/h4h-cell-cell-interactions/"

source "$HOME/miniforge3/bin/activate" "standard_env"

gtimeout 10 Rscript "$work_dir/scripts/test_timer.R" \
    --sleep_time 100

echo $?