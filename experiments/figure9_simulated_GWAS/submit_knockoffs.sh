#!/bin/bash
#
# Class: SLURM dispatcher
#
# Generate knockoffs
# 
# Author: Matteo Sesia
# Date:   06/17/2019

# Parameters
POP_LIST=("EUR" "EAS" "AMR" "SAS" "AFR")
#POP_LIST=("EUR" "EAS")
DENS_LIST=(0.1 0.05 0.01) # 0.1, 0.05
RES_LIST=("res_1" "res_2" "res_3" "res_4" "res_5" "res_6" "res_7")

# Slurm parameters
PART=candes,hns,stat,normal  # Partition names
MEMO=40G                     # Memory required (40G)
TIME=00-02:00:00              # Time required ()
CORE=1                       # Cores required ()

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" -n 1 -c "$CORE" -p "$PART" --time="$TIME

# Create directory for log files
LOGS=logs
mkdir -p $LOGS

OUT_DIR="/oak/stanford/groups/candes/matteo/invariant_knockoffs/knockoffs"
for POP in "${POP_LIST[@]}"; do    
  for DENS in "${DENS_LIST[@]}"; do    
    for RES in "${RES_LIST[@]}"; do    
      OUT_FILE=$OUT_DIR"/"$POP"_typed_"$DENS"_"$RES"_X.rds"
      if [[ ! -f $OUT_FILE ]]; then
        # Script to be run
        SCRIPT="knockoffs.sh $POP $DENS $RES"
        # Define job name
        JOBN="inv_kn_"$POP"_"$DENS"_"$RES"_s"$SEED
        OUTF=$LOGS"/"$JOBN".out"
        ERRF=$LOGS"/"$JOBN".err"
        # Assemble slurm order for this job
        ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT
        # Print order
        echo $ORD
        # Submit order
        #$ORD
        #./$SCRIPT
      fi
    done
  done
done
