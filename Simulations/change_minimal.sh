#!/bin/bash
# List of directories and corresponding values for EECC and SSCC
dirs=("0.15" "0.115" "0.082" "0.07" "0.06" "0.052" "0.05" "0.042")
# Loop through directories
for i in "${!dirs[@]}"; do
    dir="${dirs[i]}"
    cp "$dir"/rg_fiber* ./rg_fiber_"$dir"_1KX5.txt
done

