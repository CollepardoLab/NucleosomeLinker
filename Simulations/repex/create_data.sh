#!/bin/bash
# Sulis (Python)
#module load GCC/13.2.0 
#OpenMPI/4.1.6 SciPy-bundle/2023.11

module spider SciPy-bundle
module load GCC/13.2.0 OpenMPI/4.1.6 SciPy-bundle/2023.11


# Configuration directories
#declare -a configs=("0.15" "0.115" "0.082" "0.07" "0.06" "0.052" "0.05" "0.042")
configs=("0.073" "0.074" "0.075" "0.076" "0.077" "0.078" "0.079" "0.080")

for config in "${configs[@]}"; do
    cp random_second_half.py "${config}/"
    cp get_frames.sh "${config}/"
    cp multiple_dump_to_data.py "${config}/"
    cd "${config}"
    python random_second_half.py
    wait
    ./get_frames.sh random_120_frames_second_half.dump
    wait
    python multiple_dump_to_data.py data.txt
    wait
    rm frame*
    cd ../
done

wait

