#!/bin/bash

DENSITY=0.01
RES="res_1"
AMPLITUDE=20
SEED=1

DENSITY=$1
RES=$2
AMPLITUDE=$3
SEED=$4

POP_LIST=("EUR" "EAS" "AMR" "SAS" "AFR" "merged")
#POP_LIST=("merged")

OUT_DIR="/oak/stanford/groups/candes/matteo/invariant_knockoffs/stats"

for POP in "${POP_LIST[@]}"; do

  OUT_FILE=$OUT_DIR"/a"$AMPLITUDE"_s"$SEED"_"$POP"_typed_"$DENS"_"$RES"_stats.txt"
  if [[ ! -f $OUT_FILE ]]; then
    echo "########################################"
    echo "Running compute_stats.R ..."
    Rscript --vanilla compute_stats.R $POP $DENSITY $RES $AMPLITUDE $SEED
    echo "########################################"
    echo ""
  fi

done
