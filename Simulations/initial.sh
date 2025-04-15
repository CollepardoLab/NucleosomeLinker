#!/bin/bash

dirs=("0.073" "0.074" "0.075" "0.076" "0.077" "0.078" "0.079" "0.080")

files=("DNA_sequence.txt" "chromatin.so" "tremd.in" "NAFlex_params.txt" "data.txt" "unmix_dumps.py")

for i in "${!dirs[@]}"; do
    dir="${dirs[i]}"
    eecc=$(python3 salt_map_E.py "$dir")
    sscc=$(python3 salt_map_A.py "$dir")
    
    mkdir -p "$dir"
    
    for file in "${files[@]}"; do
        cp "$file" "$dir/"
    done
    
    sed -i "s/EECC/$eecc/" "$dir/tremd.in"
    sed -i "s/SSCC/$sscc/" "$dir/tremd.in"
done

echo "Files copied and tremd.in updated successfully."

