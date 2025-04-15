#!/bin/bash
#
# make_sequence.sh
#
# Usage: 
#   1) Make it executable: chmod +x make_sequence.sh
#   2) Run it: ./make_sequence.sh
#
# The script will prompt for start and end IDs, then create a file named DNA_sequence.txt
# in the current directory with lines of the form:
#   # <count_of_lines>
#   <start_ID> XX
#   <start_ID+1> XX
#   ...
#   <end_ID> XX
#

# Ask user for start/end
read -p "Enter the start ID: " START
read -p "Enter the end ID: " END

# Calculate how many lines we'll have
COUNT=$(( END - START + 1 ))

# Name of the output file
OUTPUT="DNA_sequence.txt"

echo "Creating $OUTPUT from $START to $END ..."

# Write the header
echo "# $COUNT" > "$OUTPUT"

# Loop from START to END, append lines
for (( i=START; i<=END; i++ )); do
  echo "$i XX" >> "$OUTPUT"
done

echo "Done! Check the file $OUTPUT."

