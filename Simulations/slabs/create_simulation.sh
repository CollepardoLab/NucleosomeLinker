#!/bin/bash

read -p "Enter simulation name: " SIMNAME
mkdir -p "$SIMNAME"
echo "Created folder: $SIMNAME"

# Ask how many salts, then read each salt name into an array
read -p "How many salts do you want? " NSALTS
declare -a SALTS
for (( i=1; i<=NSALTS; i++ )); do
  read -p "Enter salt #$i (e.g. 0.073): " S
  SALTS+=( "$S" )
done

read -p "How many data files (N) per salt? " N

NFOLDERS= 1

declare -a FOLDERS

for (( i=1; i<=NFOLDERS; i++ )); do
  echo
  echo "=== Data folder name==="
  read -p "Folder name (e.g. 25, Control, etc): " fname

  FOLDERS+=( "$fname" )
done

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

  mapfile -t shuffled < <(printf "%s\n" "${array[@]}" | shuf)

  local selected=("${shuffled[@]:0:count}")
  local remain=("${shuffled[@]:count}")

  data_arrays["$key"]=$(printf "%s\n" "${remain[@]}")

  printf "%s\n" "${selected[@]}"
}

for salt in "${SALTS[@]}"; do
  for (( i=0; i<NFOLDERS; i++ )); do
    init_data_array "${FOLDERS[$i]}" "$salt"
  done
done


bash create_lammps_code.sh


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

files_to_copy_parent=(
  "chromatin.so"
  "NAFlex_params.txt"
  "density.py"
  "submit_all_create.sh"
  "phase_diagram.py"
  "run_slab.sh"
)

for f in "${files_to_copy_parent[@]}"; do
  cp "$f" "$SIMNAME/"
done

rndm=$RANDOM
XXX=192673  

for salt in "${SALTS[@]}"; do
  echo
  echo "Processing salt: $salt"

  chosen=()
  total_chosen=0

  for (( i=0; i<NFOLDERS; i++ )); do
    folder="${FOLDERS[$i]}"
    count=$N

    if (( count > 0 )); then
      mapfile -t picks < <(draw_random_files "$folder" "$salt" "$count")
      chosen+=( "${picks[@]}" )
      (( total_chosen += count ))
    fi
  done

  mapfile -t chosen < <(printf "%s\n" "${chosen[@]}" | shuf)

  echo "  -> Picked $total_chosen files for salt '$salt' (ideal was $N)."


  dir="$SIMNAME/$salt"
  mkdir -p "$dir"
  mkdir -p "$dir/data"

  for (( i=0; i<total_chosen; i++ )); do
    old_path="${chosen[$i]}"
    new_id=$(( i + 1 ))
    new_file="new_data_$new_id.txt"
    cp "$old_path" "$dir/data/$new_file"

  done

  for f in "${files_to_copy[@]}"; do
    cp "$f" "$dir/"
  done

  eecc=$(python3 salt_map_E.py "$salt")
  sscc=$(python3 salt_map_A.py "$salt")
  echo "  -> eecc=$eecc, sscc=$sscc"

  sed -i "s/EECC/$eecc/" "$dir/create_slab.in"
  sed -i "s/EECC/$eecc/" "$dir/run_slab_mixture.in"

  sed -i "s/SSCC/$sscc/" "$dir/create_slab.in"
  sed -i "s/SSCC/$sscc/" "$dir/run_slab_mixture.in"

  sed -i "s/XXXX/$rndm/" "$dir/run_slab_mixture.in"
done

echo
echo "Done! Each salt in '$SIMNAME' now has up to $N unique data files, plus snippet files with sed replacements."

