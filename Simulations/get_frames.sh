#!/bin/bash

if [ -z "$1" ]; then
  echo "Usage: $0 <dump_file>"
  exit 1
fi

input_file="$1"
num_atoms=$(head -n 10 "$input_file" | grep -A1 "NUMBER OF ATOMS" | tail -1)

frame_count=0

while IFS= read -r line; do
  if [[ "$line" == "ITEM: TIMESTEP" ]]; then
    if [ $frame_count -ne 0 ]; then
      exec 3>&-
    fi
    frame_count=$((frame_count + 1))
    output_file="frame_${frame_count}.dump"
    exec 3>"$output_file"
  fi
  echo "$line" >&3
done < "$input_file"

exec 3>&-

echo "Frames have been extracted to individual files (e.g., frame_1.dump, frame_2.dump, ...)"


