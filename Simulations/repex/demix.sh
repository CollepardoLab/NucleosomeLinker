#!/bin/bash
module spider SciPy-bundle
module load GCC/13.2.0 OpenMPI/4.1.6 SciPy-bundle/2023.11
# List of directories and corresponding values for EECC and SSCC
#dirs=("0.15" "0.115" "0.082" "0.07" "0.06" "0.052" "0.05" "0.042")
dirs=("0.073" "0.074" "0.075" "0.076" "0.077" "0.078" "0.079" "0.080")
# Loop through directories
for i in "${!dirs[@]}"; do
    dir="${dirs[i]}"
    cd "$dir"/ && python unmix_dumps.py
    cd ../
done


