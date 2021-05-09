#!/bin/bash
#
# Class: SLURM dispatcher
#
# Compute test statistics
#
# Author: Matteo Sesia
# Date:   02/19/2021

# Parameters
#DENS_LIST=(0.1 0.05 0.01)
DENS_LIST=(0.01)
DENS_LIST=(0.05 0.1)
#RES_LIST=("res_1" "res_2" "res_3" "res_4" "res_5" "res_7")
RES_LIST=("res_1" "res_5") # 
AMP_LIST=(10 20 30 40 50 60 70 80 90 100)
#AMP_LIST=(10 20 50 70 100)
SEED_LIST=$(seq 1 100)

# Slurm parameters
PART=candes,stat,normal,owners  # Partition names
MEMO=10G                     # Memory required (10G)
TIME=00-06:00:00              # Time required ()
CORE=1                       # Cores required ()

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" -n 1 -c "$CORE" -p "$PART" --time="$TIME

# Create directory for log files
LOGS=logs
mkdir -p $LOGS

OUT_DIR="/oak/stanford/groups/candes/matteo/invariant_knockoffs/stats"
mkdir -p $OUT_DIR

ERASE=0
if [[ $ERASE == 1 ]]; then
  rm -f $OUT_DIR"/*.txt"
fi

for SEED in $SEED_LIST; do
  for AMPLITUDE in "${AMP_LIST[@]}"; do
    for DENS in "${DENS_LIST[@]}"; do
      for RES in "${RES_LIST[@]}"; do
        OUT_FILE=$OUT_DIR"/a"$AMPLITUDE"_s"$SEED"_merged_typed_"$DENS"_"$RES"_stats.txt"
        if [[ ! -f $OUT_FILE ]]; then
          # Script to be run
          SCRIPT="analysis.sh $DENS $RES $AMPLITUDE $SEED"
          # Define job name
          JOBN="inv_an_"$DENS"_"$RES"_s"$SEED"_a"$AMPLITUDE
          OUTF=$LOGS"/"$JOBN".out"
          ERRF=$LOGS"/"$JOBN".err"
          # Assemble slurm order for this job
          ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT
          # Print order
          echo $ORD
          # Submit order
          #$ORD
          #./$SCRIPT
          #fi
        fi
      done
    done
  done
done
