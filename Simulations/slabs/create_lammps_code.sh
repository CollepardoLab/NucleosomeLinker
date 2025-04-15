#!/bin/bash
#
# create_datafiles_in.sh
#
# Purpose:
#   Produce a SINGLE file data_files.in with:
#     1) "read_data new_data_1.txt extra/atom/types 2 group temp"
#     2) An automatically computed "change_box all x final XLO XHI y final YLO YHI z final ZLO ZHI"
#        where we expand by ±500 in x,y and ±1000 in z.
#     3) N lines of "read_data new_data_i.txt add append shift x y z"
#        using a 5x5 (x,y) grid repeated in z steps of +1000 every 25 files.
#     4) "delete_atoms group temp", "write_dump", "write_data", "write_restart"
#     5) region r1, r2 => top/bottom 100 slices in z
#
# Usage:
#   chmod +x create_datafiles_in.sh
#   ./create_datafiles_in.sh
#   Then "include data_files.in" in your LAMMPS input.
#

###############################################################################
# 1) Ask user for N
###############################################################################

read -p "How many data files (N)? " N
OUTFILE="data_files.in"
echo "Will generate $OUTFILE for N=$N"

###############################################################################
# 2) Define 5 x-values, 5 y-values => 25 combos, each block = +1000 in z
###############################################################################

xArray=(-2500 -1500 -500 500 1500)
yArray=(-2000 -1000 0 1000 2000)

zStart=-1500
zStep=1000

numXY=25  # 5 x * 5 y

###############################################################################
# 3) We need to compute all shifts (x_i, y_i, z_i) for i=1..N to find min,max.
#    We'll store them in arrays so we can do a 2-pass approach:
#      - pass 1: compute & track min/max
#      - pass 2: output lines
###############################################################################

declare -a SHIFT_LINES  # each element will be "i  xVal  yVal  zVal"

xMin=99999999
xMax=-99999999
yMin=99999999
yMax=-99999999
zMin=99999999
zMax=-99999999

for ((i=1; i<=N; i++)); do
  # block index for z
  blockIndex=$(( (i-1) / numXY ))
  offset=$(( (i-1) % numXY ))

  zVal=$(( zStart + blockIndex*zStep ))

  # figure out x,y from offset
  xIndex=$(( offset / 5 ))
  yIndex=$(( offset % 5 ))

  xVal=${xArray[$xIndex]}
  yVal=${yArray[$yIndex]}

  # store "i xVal yVal zVal"
  SHIFT_LINES+=( "$i $xVal $yVal $zVal" )

  # track min/max
  (( xVal < xMin )) && xMin=$xVal
  (( xVal > xMax )) && xMax=$xVal

  (( yVal < yMin )) && yMin=$yVal
  (( yVal > yMax )) && yMax=$yVal

  (( zVal < zMin )) && zMin=$zVal
  (( zVal > zMax )) && zMax=$zVal
done

###############################################################################
# 4) Expand bounding box:
#    x in [xMin-500, xMax+500]
#    y in [yMin-500, yMax+500]
#    z in [zMin-1000, zMax+1000]
###############################################################################

XLO=$(( xMin - 500 ))
XHI=$(( xMax + 500 ))
YLO=$(( yMin - 500 ))
YHI=$(( yMax + 500 ))
ZLO=$(( zMin - 1000 ))
ZHI=$(( zMax + 1000 ))

echo "Auto bounding box:"
echo "  x: $XLO to $XHI"
echo "  y: $YLO to $YHI"
echo "  z: $ZLO to $ZHI"

###############################################################################
# 5) region r1 => top 100 slice in z => [ZHI-100, ZHI]
#    region r2 => bottom 100 slice => [ZLO, ZLO+100]
###############################################################################

z1Lo=$(( ZHI - 150 ))
z1Hi=$(( ZHI - 50 ))

z2Lo=$(( ZLO + 50 ))
z2Hi=$(( ZLO + 150 ))

x0Lo=$(( XLO + 50 ))
x0Hi=$(( XHI - 50 ))

y0Lo=$(( YLO + 50 ))
y0Hi=$(( YHI - 50 ))

###############################################################################
# 6) Write out data_files.in with all desired lines
###############################################################################

echo "Writing $OUTFILE ..."

{
  # 6a) The initial lines:
  echo "read_data data/new_data_1.txt extra/atom/types 2 group temp"
  echo "change_box all x final $XLO $XHI y final $YLO $YHI z final $ZLO $ZHI"
  echo

  # 6b) Our repeated lines for i=1..N
  for entry in "${SHIFT_LINES[@]}"; do
    # each entry is "i xVal yVal zVal"
    set -- $entry
    i="$1"
    xVal="$2"
    yVal="$3"
    zVal="$4"
    echo "read_data data/new_data_${i}.txt add append shift $xVal $yVal $zVal"
  done

  echo
  # 6c) The final lines:
  echo "delete_atoms group temp"
  echo "write_dump all atom initial_state.dump"
  echo "write_data initial.config"
  echo "write_restart initial.restart"
  echo
  echo "region r1 block $x0Lo $x0Hi $y0Lo $y0Hi $z1Lo $z1Hi"
  echo "region r2 block $x0Lo $x0Hi $y0Lo $y0Hi $z2Lo $z2Hi"
} > "$OUTFILE"

echo "Done! You can 'include $OUTFILE' in your LAMMPS input."

