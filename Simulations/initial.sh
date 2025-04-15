#!/bin/bash

# List of directories and corresponding values for EECC and SSCC
#dirs=("0.15" "0.115" "0.082" "0.07" "0.06" "0.052" "0.05" "0.042")
dirs=("0.073" "0.074" "0.075" "0.076" "0.077" "0.078" "0.079" "0.080")
#for 1KX5 
#eecc_values=(0.4 0.375 0.35 0.3 0.25 0.2 0.1 0.01)
#sscc_values=(24 24.5 25 26 27 28 30 34)
#for H2AZ
# Files to copy
files=("DNA_sequence.txt" "chromatin.so" "tremd.in" "NAFlex_params.txt" "data.txt" "unmix_dumps.py")

# Loop through directories
for i in "${!dirs[@]}"; do
    dir="${dirs[i]}"
    eecc=$(python3 salt_map_E.py "$dir")
    sscc=$(python3 salt_map_A.py "$dir")
    
    # Create the directory if it doesn't exist
    mkdir -p "$dir"
    
    # Copy the files
    for file in "${files[@]}"; do
        cp "$file" "$dir/"
    done
    
    # Modify tremd.in
    sed -i "s/EECC/$eecc/" "$dir/tremd.in"
    sed -i "s/SSCC/$sscc/" "$dir/tremd.in"
done

echo "Files copied and tremd.in updated successfully."

