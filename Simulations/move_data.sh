#!/bin/bash

# This script should be run from Replica/, i.e. the folder where you see:
# fibers/  move_data.sh  slabs/

SOURCE="fibers"
TARGET="slabs/data"

# Folders to process under fibers/
folders=("CENPB-172")

# 1) Step into the fibers directory
cd "$SOURCE" || { echo "Could not enter $SOURCE"; exit 1; }

# 2) Loop over the top-level folders
for top in "${folders[@]}"; do
  echo "Processing $top ..."

  # Enter each top-level folder (e.g. fibers/CENPB-172)
  cd "$top" || { echo "Could not enter $top"; exit 1; }

  # 3) Look for salt directories named "0.*"
  for salt_dir in 0.*; do
    # Ensure it's actually a directory
    [ -d "$salt_dir" ] || continue

    echo "  Found salt directory: $salt_dir"

    # Create the matching directory under slabs/data/<top>/<salt_dir>
    # We are currently at: Replica/fibers/<top>/<salt_dir>
    # So, "../../$TARGET" goes back up to Replica/, then into slabs/data/
    mkdir -p "../../$TARGET/$top/$salt_dir"

    # 4) Move any new_data_* files (ignore errors if none are found)
    cp "$salt_dir"/new_data_* "../../$TARGET/$top/$salt_dir/" 2>/dev/null
  done

  # Go back up one level to "fibers/"
  cd ..
done

# Finally, go back to "Replica/"
cd ..
echo "All done!"

