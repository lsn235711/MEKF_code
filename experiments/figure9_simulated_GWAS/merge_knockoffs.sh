#!/bin/bash

DENSITY=0.01
RES="res_1"

DENSITY=$1
RES=$2

MERGE_KNOCKOFFS=1

if [[ $MERGE_KNOCKOFFS == 1 ]]; then
  echo "########################################"
  echo "Running merge_knockoffs.R ..."
  Rscript --vanilla merge_knockoffs.R $DENSITY $RES
  echo "########################################"
  echo ""
fi
