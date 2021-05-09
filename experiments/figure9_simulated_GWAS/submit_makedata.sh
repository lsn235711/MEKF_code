#!/bin/bash
#
# Class: SLURM dispatcher
#
# Generate knockoffs
# 
# Author: Matteo Sesia
# Date:   06/17/2019

# Parameters
#K_LIST=("100") # "200") #
#CHR_LIST=$(seq 1 22)

# Slurm parameters
PART=candes,hns,stat,normal,pilanci  # Partition names
MEMO=20G                     # Memory required (20G)
TIME=00-12:00:00             # Time required ()
CORE=1                       # Cores required ()

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" -n 1 -c "$CORE" -p "$PART" --time="$TIME

# Create directory for log files
LOGS=logs
mkdir -p $LOGS

OUT_DIR="/oak/stanford/groups/candes/matteo/invariant_knockoffs/"
#    for CHR in $CHR_LIST; do
#      for RES in "${RES_LIST[@]}"; do    
#OUT_FILE=$OUT_DIR"/ukb_gen_chr"$CHR"_ibd"$IBD"_res"$RES".bed"
#if [[ ! -f $OUT_FILE ]]; then
# Script to be run
SCRIPT="make_data.sh"
# Define job name
JOBN="inv_dat"
OUTF=$LOGS"/"$JOBN".out"
ERRF=$LOGS"/"$JOBN".err"
# Assemble slurm order for this job
ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT
# Print order
echo $ORD
# Submit order
$ORD
#./$SCRIPT
#fi
#      done
#    done
