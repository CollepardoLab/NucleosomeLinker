#!/bin/bash

# This script should be run from Replica/, i.e. the folder where you see:
# repex/  move_data.sh  slabs/

SOURCE="repex"
TARGET="slabs/data"

# Folders to process under fibers/
folders=("25")

cd "$SOURCE" || { echo "Could not enter $SOURCE"; exit 1; }

for top in "${folders[@]}"; do
  echo "Processing $top ..."

  cd "$top" || { echo "Could not enter $top"; exit 1; }

  for salt_dir in 0.*; do
  
    [ -d "$salt_dir" ] || continue

    echo "  Found salt directory: $salt_dir"

    # Create the matching directory under slabs/data/<top>/<salt_dir>
    # We are currently at: repex/<top>/<salt_dir>
    # So, "../../$TARGET" goes back up to Parent/, then into slabs/data/
    mkdir -p "../../$TARGET/$top/$salt_dir"

    cp "$salt_dir"/new_data_* "../../$TARGET/$top/$salt_dir/" 2>/dev/null
  done

  cd ..
done

cd ..
echo "All done!"

