#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

density <- 0.01
res <- "res_1"

density <- as.numeric(args[1])
res <- as.character(args[2])

suppressMessages(library(tidyverse))
library(bigstatsr)

idir <- "/oak/stanford/groups/candes/matteo/invariant_knockoffs"

pop.list <- c("EUR", "EAS", "AMR", "SAS", "AFR")

## Load knockoff legend info
ifile <- sprintf("%s/knockoffs/%s_typed_%s_%s_legend.txt", idir, pop.list[1], density, res)
legend <- read_delim(ifile, delim=" ", col_types=cols())

## Load knockoff sample info
sample <- do.call("rbind", lapply(pop.list, function(pop) {
  ifile <- sprintf("%s/knockoffs/%s_typed_%s_%s_sample.txt", idir, pop, density, res)
  sample.pop <- read_delim(ifile, delim=" ", col_types=cols())
}))

## Initialize storage for merged knockoffs
ofile <- sprintf("%s/knockoffs/merged_typed_%s_%s_X", idir, density, res)
fbm <- FBM(nrow(sample), nrow(legend), type = c("unsigned short"), backingfile = ofile)

for(pop in pop.list) {
  cat(sprintf("Copying data for %s population...\n", pop))
  ## Load knockoff-augmented data
  ifile <- sprintf("%s/knockoffs/%s_typed_%s_%s_X.rds", idir, pop, density, res)
  fbm.pop <- bigstatsr::big_attach(ifile)
  ## Copy the data
  idx.rows <- which(sample$population==pop)
  fbm[idx.rows,] <- fbm.pop[,]
}
cat(sprintf("Saving merged data...\n", pop))
fbm$save()
cat(sprintf("Output file: %s\n", ofile))

ofile <- sprintf("%s/knockoffs/merged_typed_%s_%s_sample.txt", idir, density, res)
sample %>% write_delim(ofile, delim=" ")

ofile <- sprintf("%s/knockoffs/merged_typed_%s_%s_legend.txt", idir, density, res)
legend %>% write_delim(ofile, delim=" ")
