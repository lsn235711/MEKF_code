suppressMessages(library(tidyverse))
library(bigstatsr)

idir <- "/oak/stanford/groups/candes/matteo/invariant_knockoffs/data"
odir <- "/oak/stanford/groups/candes/matteo/invariant_knockoffs/data"

## Load sample file
ifile <- sprintf("%s/full_sample.txt", idir)
sample <- read_delim(ifile, delim=" ", col_types=cols())

## Load legend file
ifile <- sprintf("%s/full_legend.txt", idir)
legend <- read_delim(ifile, delim=" ", col_types=cols())

## Load genotypes
ifile <- sprintf("%s/full_X.rds", idir)
X.fbm <- bigstatsr::big_attach(ifile)

## Load haplotypes
ifile <- sprintf("%s/full_H.rds", idir)
H.fbm <- bigstatsr::big_attach(ifile)

###############
## Type data ##
###############

source("utils.R")

## Compute typing indices
density.list <- c(0.1, 0.05, 0.01)

cat(sprintf("Computing typing indices... "))
idx.typed.list <- compute_typing_indices(legend, density.list)
cat("done.\n")

cat(sprintf("Saving data... "))

for(k in 1:length(density.list)) {
  density <- density.list[k]
  idx.typed <- idx.typed.list[[k]]
  
  ## Save typed haplotype data
  ofile <- sprintf("%s/typed_%s_H", odir, density)
  H.fbm.typed <- big_copy(H.fbm, ind.col=idx.typed, backingfile=ofile)
  H.fbm.typed$save()

  ## Save typed genotype data
  ofile <- sprintf("%s/typed_%s_X", odir, density)
  X.fbm.typed <- big_copy(X.fbm, ind.col=idx.typed, backingfile=ofile)
  X.fbm.typed$save()
  
  ## Save legend file
  ofile <- sprintf("%s/typed_%s_legend.txt", odir, density)
  legend[idx.typed,] %>% write_delim(ofile, delim=" ")
  
  ## Save sample file
  ofile <- sprintf("%s/typed_%s_sample.txt", odir, density)
  sample %>% write_delim(ofile, delim=" ")
  
}

cat("done.\n")
cat(sprintf("Data are stored in:\n  %s\n", odir))
