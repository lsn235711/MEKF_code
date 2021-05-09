#!/bin/bash

# Parameters
PHENO_LIST=("platelet" "height" "bmi" "cvd" "diabetes" "hypothyroidism" "respiratory" "sbp")
RES_LIST=("res6" "res5" "res4" "res3" "res2" "res1" "res0")
R_LIST=(5)
Q_LIST=(0.1 0.2)
NREPS=100

# Slurm parameters
MEMO=2G                             # Memory required (2GB)
TIME=02:00:00                       # Time required (30m)

# Assemble order prefix
ORDP="sbatch --mem="$MEMO" --nodes=1 --ntasks=1 --cpus-per-task=1 --time="$TIME

# Create directory for log files
LOGS="logs"
mkdir -p $LOGS

OUT_DIR="results"
mkdir -p $OUT_DIR

# Loop over configurations and chromosomes
for PHENO in "${PHENO_LIST[@]}"; do
  for RES in "${RES_LIST[@]}"; do
    for Q in "${Q_LIST[@]}"; do
      for R in "${R_LIST[@]}"; do

        JOBN=$PHENO"_"$RES"_q"$Q"_r"$R

        OUT_FILE=$OUT_DIR"/discoveries_"$JOBN".txt"

        RUN=1
        if [[ -f $OUT_FILE ]]; then
          RUN=1
          #echo "Found results in "$OUT_FILE
        fi

        if [[ $RUN == 1 ]]; then
          # Script to be run
          SCRIPT="compare.sh $PHENO $RES $Q $R $NREPS"
          # Define job name for this chromosome
          OUTF=$LOGS"/"$JOBN".out"
          ERRF=$LOGS"/"$JOBN".err"
          # Assemble slurm order for this job
          ORD=$ORDP" -J "$JOBN" -o "$OUTF" -e "$ERRF" "$SCRIPT
          # Print order
          echo $SCRIPT
          # Submit order
          $ORD
          # Run command now
          #./$SCRIPT
        fi
      done
    done
  done
done
