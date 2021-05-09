#!/bin/bash

mkdir -p /oak/stanford/groups/candes/matteo/invariant_knockoffs/
mkdir -p /oak/stanford/groups/candes/matteo/invariant_knockoffs/data/
mkdir -p /oak/stanford/groups/candes/matteo/invariant_knockoffs/hmm/
mkdir -p /oak/stanford/groups/candes/matteo/invariant_knockoffs/knockoffs/
mkdir -p /oak/stanford/groups/candes/matteo/invariant_knockoffs/stats/
mkdir -p $PI_SCRATCH"/matteo"
mkdir -p $PI_SCRATCH"/matteo/invariant_knockoffs"

N=50000
#N=1000

ERASE=0

GENERATE_DATA=0
PARTITION_GENOME=0
TYPE_DATA=0
GENERATE_Y=0

if [[ $ERASE == 1 ]]; then
  DATA_DIR="/scratch/groups/candes/matteo/invariant_knockoffs"
  rm -f $DATA_DIR/data/*.bk
  rm -f $DATA_DIR/data/*.rds
  rm -f $DATA_DIR/hmm/*.bk
  rm -f $DATA_DIR/hmm/*.rds
fi

if [[ $GENERATE_DATA == 1 ]]; then
  echo "########################################"
  echo "Running generate_data.R ..."
  Rscript --vanilla generate_data.R $N
  echo "########################################"
  echo ""
fi

if [[ $PARTITION_GENOME == 1 ]]; then
  echo "########################################"
  echo "Running partition_genome.R ..."
  Rscript --vanilla partition_genome.R
  echo "########################################"
  echo ""
fi

if [[ $TYPE_DATA == 1 ]]; then
  echo "########################################"
  echo "Running type_data.R ..."
  Rscript --vanilla type_data.R
  echo "########################################"
  echo ""
fi

if [[ $GENERATE_Y == 1 ]]; then
  echo "########################################"
  echo "Running generate_phenotypes.R ..."
  Rscript --vanilla generate_phenotypes.R
  echo "########################################"
  echo ""
fi
