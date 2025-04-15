#!/bin/bash

# List of directories and corresponding values for EECC and SSCC
dirs=("0.15" "0.115" "0.082" "0.07" "0.06" "0.052" "0.05" "0.042")
#for 1KX5 
#eecc_values=(0.4 0.375 0.35 0.3 0.25 0.2 0.1 0.01)
#sscc_values=(24 24.5 25 26 27 28 30 34)
# For H2AZ
eecc_values=(0.35 0.325 0.3 0.25 0.2 0.15 0.05 0.005)
sscc_values=(23 23.5 24 25 26 27 29 33)
rndm=$RANDOM
# Files to copy
files=("DNA_sequence.txt" "chromatin.so" "in.create_slab2" "NAFlex_params.txt" "run_slab_mixture.in")

# Loop through directories
for i in "${!dirs[@]}"; do
    dir="${dirs[i]}"
    eecc="${eecc_values[i]}"
    sscc="${sscc_values[i]}"
    
    # Create the directory if it doesn't exist
    mkdir -p "$dir"
    cp ../../../fibers/H2AZ_fiber_for_slab/"$dir"/3new_data_*.txt "$dir/"
    # Copy the files
    for file in "${files[@]}"; do
        cp "$file" "$dir/"
    done
    
    # Modify tremd.in
    sed -i "s/EECC/$eecc/" "$dir/in.create_slab2"
    sed -i "s/EECC/$eecc/" "$dir/run_slab_mixture.in"
    sed -i "s/SSCC/$sscc/" "$dir/in.create_slab2"
    sed -i "s/SSCC/$sscc/" "$dir/run_slab_mixture.in"
    sed -i "s/XXXX/$rndm/" "$dir/run_slab_mixture.in"
done


