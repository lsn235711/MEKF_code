#!/bin/bash

POP="EUR"
DENSITY=0.01
RES="res_1"

POP=$1
DENSITY=$2
RES=$3

mkdir -p /oak/stanford/groups/candes/matteo/invariant_knockoffs/knockoffs/

ERASE=0

GENERATE_KNOCKOFFS=1

if [[ $ERASE == 1 ]]; then
  rm -f /oak/stanford/groups/candes/matteo/invariant_knockoffs/knockoffs/*.bk
  rm -f /oak/stanford/groups/candes/matteo/invariant_knockoffs/knockoffs/*.rds
fi

if [[ $GENERATE_KNOCKOFFS == 1 ]]; then
  echo "########################################"
  echo "Running generate_knockoffs.R ..."
  Rscript --vanilla generate_knockoffs.R $POP $DENSITY $RES
  echo "########################################"
  echo ""
fi
