#!/bin/bash

# Check if a file is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <dump_file>"
  exit 1
fi

# Input file and number of atoms
input_file="$1"
num_atoms=$(head -n 10 "$input_file" | grep -A1 "NUMBER OF ATOMS" | tail -1)

# Frame counter
frame_count=0

# Read the file line by line
while IFS= read -r line; do
  # Detect the start of a new frame
  if [[ "$line" == "ITEM: TIMESTEP" ]]; then
    # Close the previous frame file if it exists
    if [ $frame_count -ne 0 ]; then
      exec 3>&-
    fi
    # Increment frame counter and open a new file for the next frame
    frame_count=$((frame_count + 1))
    output_file="frame_${frame_count}.dump"
    exec 3>"$output_file"
  fi
  # Write line to the current frame file
  echo "$line" >&3
done < "$input_file"

# Close the last frame file
exec 3>&-

echo "Frames have been extracted to individual files (e.g., frame_1.dump, frame_2.dump, ...)"


