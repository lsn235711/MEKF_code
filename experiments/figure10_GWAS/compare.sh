#!/bin/bash

PHENOTYPE=$1
RESOLUTION=$2
Q=$3
R=$4
NREPS=$5

ml gcc/8.3.0
ml r/4.0.0

Rscript --vanilla compare_discoveries.R $PHENOTYPE $RESOLUTION $Q $R $NREPS
