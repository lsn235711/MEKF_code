#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

pop <- "EUR"
density <- 0.01
res <- "res_1"

pop <- as.character(args[1])
density <- as.numeric(args[2])
res <- as.character(args[3])

seed <- 2021

suppressMessages(library(tidyverse))
library(SNPknock)
library(bigstatsr)

idir <- "/oak/stanford/groups/candes/matteo/invariant_knockoffs"
odir <- "/oak/stanford/groups/candes/matteo/invariant_knockoffs/knockoffs"

source("utils.R")

## Load sample file
ifile <- sprintf("%s/data/full_sample.txt", idir)
sample <- read_delim(ifile, delim=" ", col_types=cols())
populations <- unique(sample$population)

## Load full legend file
ifile <- sprintf("%s/data/full_legend.txt", idir)
legend.full <- read_delim(ifile, delim=" ", col_types=cols())

## Load partition file
ifile <- sprintf("%s/data/legend_partitions.txt", idir)
partitions <- read_delim(ifile, delim=" ", col_types=cols())

## Load legend file
ifile <- sprintf("%s/data/typed_%s_legend.txt", idir, density)
legend <- read_delim(ifile, delim=" ", col_types=cols())

## Groups at this resolution
idx.typed <- which(legend.full$id %in% legend$id)
partition <- partitions[idx.typed,]
partition$group <- partition[[res]]
partition <- partition %>% select(id, position, map, group)
groups <- cumsum(c(1,diff(partition$group)>0))
storage.mode(groups) <- "integer"

## Load haplotypes
ifile <- sprintf("%s/data/typed_%s_H.rds", idir, density)
H <- bigstatsr::big_attach(ifile)

## Specify population
idx.pop <- which(sample$population==pop)

## Load HMM
ifile <- sprintf("%s/hmm/%s_H.rds", idir, pop)
fbm.ref <- bigstatsr::big_attach(ifile)
hmm <- initialize_hmm_single(fbm.ref, legend, legend.full)

## Number of individuals in this population
n <- length(idx.pop)/2
## Number of variants
p <- ncol(H)

## Randomly swap genotypes and knockoffs
cat("Generating random swaps... ")
set.seed(seed)
n.groups <- max(groups)
group.swap <- sort(sample(n.groups, round(n.groups/2)))
idx.swap <- which(groups %in% group.swap)
cat("done\n")

## Initialize knockoff container (too slow)
#ofile <- sprintf("%s/%s_typed_%s_%s_X", odir, pop, density, res)
#X.all <- FBM(n, 2*p, type="unsigned short", backingfile=ofile)

X.all.buf <- matrix(0, n, 2*p)
  
## Generate knockoffs
cat("Generating knockoffs... ")
pb <- txtProgressBar(min = 0, max = n, style = 3, file=stderr())
x.all <- rep(0, 2*p)
for(i in 1:n) {
  i.idx <- idx.pop[c(2*i-1,2*i)]
  hk <- knockoffHaplotypes(H[i.idx,,drop=FALSE], hmm$r, hmm$alpha, hmm$theta, groups=groups, seed=i)
  xk <- hk[1,] + hk[2,]
  x <- H[i.idx[1],] + H[i.idx[2],]
  if(length(idx.swap)>0) {
    x.all[-c(idx.swap,idx.swap+p)] <- c(x[-idx.swap], xk[-idx.swap])
    x.all[c(idx.swap,idx.swap+p)] <- c(xk[idx.swap], x[idx.swap])
  } else {
    x.all <- c(x,xk)
  }
  X.all.buf[i,] <- x.all
  setTxtProgressBar(pb, i)
}
close(pb)
cat("done\n")

## Save knockoff-augmented genotypes
cat("Saving knockoffs... ")
ofile <- sprintf("%s/%s_typed_%s_%s_X", odir, pop, density, res)
fbm <- bigstatsr::as_FBM(X.all.buf, type = "unsigned short", backingfile=ofile)
fbm$save()
cat("done\n")

## Save legend
legend.1 <- legend %>% mutate(knockoff=FALSE, group=groups) %>%
  select(id, position, knockoff, group, everything())
legend.2 <- legend %>% mutate(knockoff=TRUE, group=groups) %>%
  select(id, position, knockoff, group, everything())
legend.1$knockoff[idx.swap] <- !legend.1$knockoff[idx.swap]
legend.2$knockoff[idx.swap] <- !legend.2$knockoff[idx.swap]
legend.out <- rbind(legend.1, legend.2)
ofile <- sprintf("%s/%s_typed_%s_%s_legend.txt", odir, pop, density, res)
legend.out %>% write_delim(ofile, delim=" ")

## Save sample information
ofile <- sprintf("%s/%s_typed_%s_%s_sample.txt", odir, pop, density, res)
sample[idx.pop[seq(1,length(idx.pop),by=2)],] %>% write_delim(ofile, delim=" ")

cat(sprintf("Knockoffs are stored in:\n  %s\n", odir))
