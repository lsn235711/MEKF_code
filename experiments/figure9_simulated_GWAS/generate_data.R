#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

n <- 100

## Number of samples to generate
n <- as.integer(args[1])

suppressMessages(library(tidyverse))
library(bigstatsr)

idir <- "/oak/stanford/groups/candes/matteo/invariant_knockoffs/data_1000gp/1000GP_Phase3_clean"
##odir <- "/oak/stanford/groups/candes/matteo/invariant_knockoffs"
odir <- "/scratch/groups/candes/matteo/invariant_knockoffs"

## Load sample file
ifile <- sprintf("%s/1000GP_Phase3.sample", idir)
sample <- read_delim(ifile, delim=" ", col_types=cols())

## Load legend file
ifile <- sprintf("%s/legend_chr22.txt", idir)
legend <- read_delim(ifile, delim=" ", col_types=cols())

## Import haplotypes from FBM
rds.file <- sprintf("%s/1000GP_Phase3_chr22.rds", idir, ifile)
fbm <- bigstatsr::big_attach(rds.file)

## Load useful functions
source("utils.R")

## Create reference panel
cat(sprintf("Creating reference panel... "))
populations <- unique(sample$GROUP)
ref1000gp <- subset_1000gp(fbm, sample, legend, populations, M=20)
cat("done.\n")

## Initialize HMM
cat(sprintf("Initializing HMM... "))
hmm.list <- initialize_hmm(ref1000gp)
cat("done.\n")

## Generate haplotype data
cat(sprintf("Generating haplotype and genotype data... \n"))

## Initialize FBMs
p <- nrow(legend)
ofile <- sprintf("%s/data/full_H", odir)
H <- FBM(2*n, p, type="unsigned short", backingfile=ofile)
ofile <- sprintf("%s/data/full_X", odir)
X <- FBM(n, p, type="unsigned short", backingfile=ofile)

## generate_genotypes <- function(hmm.list, legend, n, return.Z=FALSE) {
##   K <- length(hmm.list)
##   populations.int <- rep(1:K, each=ceiling(2*n/K))[1:(2*n)]
##   pop.n <- sapply(1:K, function(k) sum(populations.int==k))
##   block.size <- 10
##   for(k in 1:K) {
##     cat(sprintf("Generating haplotype and genotype data for population %d of %d... \n", k, K))
##     hmm <- hmm.list[[k]]
    
##     pop.idx <- which(populations.int==k)

##     idx.row.list <- split(pop.idx, ceiling(seq_along(pop.idx)/block.size))
##     for(b in 1:length(idx.row.list)) {
##       cat(sprintf("Processing block %d of %d...\n", b, length(idx.row.list)))
##       idx.row <- idx.row.list[[b]]
##       idx.row.x <- round(idx.row[seq(2,length(idx.row),by=2)]/2)
##       indices.h1 <- seq(1,length(idx.row),by=2)
##       indices.h2 <- seq(2,length(idx.row),by=2)
      
##       h <- generate_haplotypes_onepop(hmm, idx.row)
      
##       cat(sprintf("Storing data for block %d of %d... ", b, length(idx.row.list)))
##       H[idx.row,] <- h
##       X[idx.row.x,] <- h[indices.h1,] + h[indices.h2,]
##       cat("done.\n")
##     }
##   }
##   out <- c()
##   out$legend <- legend
##   out$sample <- tibble(id=rep(1:n, each=2), hap=rep(c(0,1), length.out=2*n),
##                        population=populations[populations.int])
##   return(out)
## }
generate_genotypes <- function(hmm.list, legend, n, return.Z=FALSE) {
  K <- length(hmm.list)
  populations.int <- rep(1:K, each=ceiling(2*n/K))[1:(2*n)]
  pop.n <- sapply(1:K, function(k) sum(populations.int==k))
  for(k in 1:K) {
    cat(sprintf("Generating haplotype and genotype data for population %d of %d... \n", k, K))
    hmm <- hmm.list[[k]]
    idx.row <- which(populations.int==k)
    generate_haplotypes_onepop(hmm, idx.row)      
  }
  out <- c()
  out$legend <- legend
  out$sample <- tibble(id=rep(1:n, each=2), hap=rep(c(0,1), length.out=2*n),
                       population=populations[populations.int])
  return(out)
}

data.full <- generate_genotypes(hmm.list, legend, n)

###############
## Save data ##
###############

cat(sprintf("Saving data... "))

## Save legend
ofile <- sprintf("%s/data/full_legend.txt", odir)
data.full$legend %>% write_delim(ofile, delim=" ")

## Save list of samples
ofile <- sprintf("%s/data/full_sample.txt", odir)
data.full$sample %>% write_delim(ofile, delim=" ")

## Save haplotypes and genotype FBM
H$save()
X$save()

cat("done.\n")
cat(sprintf("Data is stored in:\n  %s/data/\n", odir))

###############
## Save HMM  ##
###############
cat(sprintf("Saving hmm... "))
for(k in 1:length(populations)) {
  pop = populations[k]
  ## Store sample information
  ofile <- sprintf("%s/hmm/%s_sample.txt", odir, pop)
  ref1000gp[[k]]$sample %>% write_delim(ofile, delim=" ")
  ## Store legend information
  ofile <- sprintf("%s/hmm/%s_legend.txt", odir, pop)
  ref1000gp[[k]]$variants %>% write_delim(ofile, delim=" ")
  ## Store reference haplotypes
  ofile <- sprintf("%s/hmm/%s_H", odir, pop)
  fbm <- bigstatsr::as_FBM(ref1000gp[[k]]$H, type = "unsigned short", backingfile=ofile)
  fbm$save()  
}
cat("done.\n")
cat(sprintf("HMM is stored in:\n  %s\n", odir))
