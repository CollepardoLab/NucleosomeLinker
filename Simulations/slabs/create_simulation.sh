#!/bin/bash
#
# create_sim.sh
#
# Purpose:
#   1) Prompt user for:
#      - Simulation name
#      - A list of salts (user-input)
#      - Number of data files N per salt
#      - Data folders & percentages (must sum to 100)
#   2) For each salt in the user-defined list:
#      - For each data folder, randomly pick (percentage*N/100) from data/<folder>/<salt>/new_data_*
#      - Copy them to SIMNAME/<salt>/data, named as new_data_1..new_data_k
#      - Copy snippet files, do sed replacements in create_slab.in & run_slab_mixture.in:
#         "s/EECC/$eecc/", "s/SSCC/$sscc/", "s/XXXX/$rndm/"
#      - eecc, sscc read from salt_map_E.py, salt_map_A.py
#
# Run from inside slabs/:  ./create_sim.sh
#

###############################################################################
# 1) Ask user for simulation name, number of salts, etc.
###############################################################################

read -p "Enter simulation name: " SIMNAME
mkdir -p "$SIMNAME"
echo "Created folder: $SIMNAME"

# Ask how many salts, then read each salt name into an array
read -p "How many salts do you have? " NSALTS
declare -a SALTS
for (( i=1; i<=NSALTS; i++ )); do
  read -p "Enter salt #$i (e.g. 0.073): " S
  SALTS+=( "$S" )
done

read -p "How many data files (N) per salt? " N

read -p "How many data folders to blend? " NFOLDERS

declare -a FOLDERS
declare -a PERCENTAGES
sum_perc=0

for (( i=1; i<=NFOLDERS; i++ )); do
  echo
  echo "=== Data folder $i of $NFOLDERS ==="
  read -p "Folder name (e.g. CENPB-172, Control): " fname
  read -p "What percentage of $N should come from '$fname'? " perc

  FOLDERS+=( "$fname" )
  PERCENTAGES+=( "$perc" )
  (( sum_perc += perc ))
done

if (( sum_perc != 100 )); then
  echo "ERROR: Percentages must sum to 100. Got $sum_perc."
  exit 1
fi

###############################################################################
# 2) We'll store each folder+salt's new_data_* in an associative array:
#    data_arrays["folder:salt"]="(list of files)"
#    draw_random_files picks and *removes* them so they aren't reused.
###############################################################################

declare -A data_arrays

init_data_array() {
  local folder="$1"
  local salt="$2"
  local path="data/$folder/$salt"

  if [[ -d "$path" ]]; then
    # Gather all new_data_* in this folder for this salt
    local all_files
    all_files=$(find "$path" -type f -name "new_data_*" 2>/dev/null)
    data_arrays["$folder:$salt"]="$all_files"
  else
    # If folder/salt subdir not found => no files
    data_arrays["$folder:$salt"]=""
  fi
}

draw_random_files() {
  local folder="$1"
  local salt="$2"
  local count="$3"

  local key="$folder:$salt"
  local list="${data_arrays["$key"]}"

  # Convert multiline string to array
  readarray -t array <<< "$list"
  local total="${#array[@]}"

  if (( total < count )); then
    echo "ERROR: $folder/$salt only has $total files, need $count."
    exit 1
  fi

  # Shuffle them
  mapfile -t shuffled < <(printf "%s\n" "${array[@]}" | shuf)

  # First 'count' are selected
  local selected=("${shuffled[@]:0:count}")
  local remain=("${shuffled[@]:count}")

  # Update the array with leftover (i.e. remove these from the pool)
  data_arrays["$key"]=$(printf "%s\n" "${remain[@]}")

  # Output the selected list so caller can consume it
  printf "%s\n" "${selected[@]}"
}

###############################################################################
# Initialize the data_arrays once for every (folder, salt).
###############################################################################
for salt in "${SALTS[@]}"; do
  for (( i=0; i<NFOLDERS; i++ )); do
    init_data_array "${FOLDERS[$i]}" "$salt"
  done
done

###############################################################################
# 3) We'll run our code snippet to prep something (per original script)
###############################################################################

bash create_lammps_code.sh

###############################################################################
# 4) Define the files that get copied:
###############################################################################

# These go inside each salt's directory
files_to_copy=(
  "data_files.in"
  "chromatin.so"
  "NAFlex_params.txt"
  "create_slab.in"
  "run_slab.sub"
  "run_slab_mixture.in"
  "run_slab.sh"
  "create_slab.sub"
  "DNA_sequence.txt"
)

# These get copied once at the top level of SIMNAME
files_to_copy_parent=(
  "chromatin.so"
  "NAFlex_params.txt"
  "density.py"
  "submit_all_create.sh"
  "phase_diagram.py"
  "run_slab.sh"
)

# Copy parent-level files now (only once)
for f in "${files_to_copy_parent[@]}"; do
  cp "$f" "$SIMNAME/"
done

###############################################################################
# 5) We'll do "s/XXXX/$rndm/" as requested, for the snippet
###############################################################################
rndm=$RANDOM
XXX=192673  # For reference, though unused in sed below

echo
echo "===== Distributing data across salts, applying snippet logic ====="

# Prepare a single CSV for the entire simulation, with one header line
csvfile="$SIMNAME/fiber_classification_per_salt.csv"
echo "Salt,FiberID,Type" > "$csvfile"

###############################################################################
# 6) For each salt, pick the fraction from each folder, rename them, copy snippet
###############################################################################
for salt in "${SALTS[@]}"; do
  echo
  echo "Processing salt: $salt"

  # We'll gather all picks for this salt in 'chosen'
  chosen=()
  total_chosen=0

  # For each folder, pick (percentage*N/100)
  for (( i=0; i<NFOLDERS; i++ )); do
    folder="${FOLDERS[$i]}"
    perc="${PERCENTAGES[$i]}"
    count=$(( (perc * N) / 100 ))

    if (( count > 0 )); then
      mapfile -t picks < <(draw_random_files "$folder" "$salt" "$count")
      chosen+=( "${picks[@]}" )
      (( total_chosen += count ))
    fi
  done

  # Shuffle the final chosen list to mix them
  mapfile -t chosen < <(printf "%s\n" "${chosen[@]}" | shuf)

  echo "  -> Picked $total_chosen mixed files for salt '$salt' (ideal was $N)."

  # Create the salt directory under SIMNAME
  dir="$SIMNAME/$salt"
  mkdir -p "$dir"
  mkdir -p "$dir/data"

  # Rename each chosen file to new_data_1..(total_chosen)
  normal_count=0
  modified_count=0
  for (( i=0; i<total_chosen; i++ )); do
    old_path="${chosen[$i]}"
    new_id=$(( i + 1 ))
    new_file="new_data_$new_id.txt"
    cp "$old_path" "$dir/data/$new_file"

    # Check the # of atoms from the file (3rd line has it in your format)
    atom_count=$(head -3 "$old_path" | tail -1 | awk '{print $1}')
    if (( atom_count > 424 )); then
      echo "$salt,$new_id,modified" >> "$csvfile"
      (( modified_count++ ))
    else
      echo "$salt,$new_id,normal" >> "$csvfile"
      (( normal_count++ ))
    fi
  done

  echo "Salt: $salt -> Normal: $normal_count, Modified: $modified_count"

  # Copy snippet files into each salt subdir
  for f in "${files_to_copy[@]}"; do
    cp "$f" "$dir/"
  done

  # Get eecc, sscc from Python
  eecc=$(python3 salt_map_E.py "$salt")
  sscc=$(python3 salt_map_A.py "$salt")
  echo "  -> eecc=$eecc, sscc=$sscc"

  # Perform sed replacements in create_slab.in & run_slab_mixture.in
  sed -i "s/EECC/$eecc/" "$dir/create_slab.in"
  sed -i "s/EECC/$eecc/" "$dir/run_slab_mixture.in"

  sed -i "s/SSCC/$sscc/" "$dir/create_slab.in"
  sed -i "s/SSCC/$sscc/" "$dir/run_slab_mixture.in"

  sed -i "s/XXXX/$rndm/" "$dir/run_slab_mixture.in"
done

echo
echo "Done! Each salt in '$SIMNAME' now has up to N unique data files, plus snippet files with sed replacements."

