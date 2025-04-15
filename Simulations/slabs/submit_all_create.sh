#!/bin/bash

# Get the current user's username
USER=$(whoami)

# Loop through all directories that look like numerical salt concentrations
for dir in */; do
    # Remove trailing slash
    dir=${dir%/}

    # Check if `create_slab.sub` exists in the directory
    if [[ -f "$dir/create_slab.sub" ]]; then
        echo "Submitting job in $dir..."
        (cd "$dir" && sbatch create_slab.sub)
    else
        echo "Skipping $dir (no create_slab.sub found)"
    fi
done

# Show submitted jobs
echo "Checking queue for user $USER..."
squeue -u "$USER"

