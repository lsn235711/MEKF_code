#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

knockoff.threshold <- function(W, fdr=0.10, offset=1) {
  if(offset>1 | offset<0) {
    stop('Input offset must be between 0 or 1')
  }
  ts = sort(c(0, abs(W)))
  ratio = sapply(ts, function(t)
    (offset + sum(W <= -t)) / max(1, sum(W >= t)))
  ok = which(ratio <= fdr)
  ifelse(length(ok) > 0, ts[ok[1]], Inf)
}

pop <- "merged"
density <- 0.1
res <- "res_1"
amplitude <- 100
seed <- 9

pop <- as.character(args[1])
density <- as.numeric(args[2])
res <- as.character(args[3])
amplitude <- as.numeric(args[4])
seed <- as.numeric(args[5])

suppressMessages(library(tidyverse))
library(bigstatsr)

idir <- "/oak/stanford/groups/candes/matteo/invariant_knockoffs"
odir <- "/oak/stanford/groups/candes/matteo/invariant_knockoffs/stats"

source("utils.R")

## Load phenotypes
ifile <- sprintf("%s/data/phenotypes.txt", idir)
phenotypes <- read_delim(ifile, delim=" ", col_types=cols())

## Load list of causal variants
ifile <- sprintf("%s/data/legend_causal.txt", idir)
legend.causal <- read_delim(ifile, delim=" ", col_types=cols())

## Load partition file
ifile <- sprintf("%s/data/legend_partitions.txt", idir)
partitions <- read_delim(ifile, delim=" ", col_types=cols())
partition <- partitions %>% select(id, position, map) %>%
  mutate(causal = (id %in% legend.causal$id))
partition$group <- partitions[[res]]
partition <- partition %>% group_by(group) %>%
  summarise(bp.min=min(position), bp.max=max(position), causal=any(causal)) %>%
  mutate(width = (bp.max-bp.min)/1e6)

## Load knockoff sample info
ifile <- sprintf("%s/knockoffs/%s_typed_%s_%s_sample.txt", idir, pop, density, res)
sample <- read_delim(ifile, delim=" ", col_types=cols())

## Load knockoff legend info
ifile <- sprintf("%s/knockoffs/%s_typed_%s_%s_legend.txt", idir, pop, density, res)
legend <- read_delim(ifile, delim=" ", col_types=cols())

## Load knockoff-augmented data
ifile <- sprintf("%s/knockoffs/%s_typed_%s_%s_X.rds", idir, pop, density, res)
fbm <- bigstatsr::big_attach(ifile)

## Compute scaling factors
scaler <- big_scale()
G.scale <- scaler(fbm)
scaling.factors <- G.scale$scale

## Fit the lasso
pheno.name <- sprintf("Y_a%s_s%s", amplitude, seed)
idx.pop <- which(phenotypes$id %in% sample$id)
y <- phenotypes[idx.pop,][[pheno.name]]
lasso.fit <- big_spLinReg(fbm, y)

## Extracting lasso coefficients
cat("- Extracting model coefficients... ")
betas <- sapply(1:10, function(k) lasso.fit[[1]][k][[1]]$beta)
beta.extracted <- get_beta(betas, method = c("mean"))
beta.nonzero <- beta.extracted[1:(length(beta.extracted))]
lasso.indices <- attr(lasso.fit, "ind.col")
beta.variants <- rep(0, ncol(fbm))
beta.variants[lasso.indices] <- beta.nonzero
model.support <- which(beta.variants!=0)
cat("OK.\n")
cat(sprintf("- Model support size = %d.\n", length(model.support)))

## Compute importance measures
cat("- Computing importance measures... ")
Lasso <- legend %>% mutate(Beta = beta.variants,
                          Scale=scaling.factors,
                          Z=Beta*Scale) %>%
  filter(Z!=0) %>%
  arrange(desc(abs(Z))) %>%
  select(id, position, knockoff, group, Beta, Scale, Z)
cat("OK.\n")

Stats <- Lasso %>%
  group_by(group, knockoff) %>%
  summarise(position.lead=position[which.max(abs(Z))], Z=sum(abs(Z))) %>%
  group_by(group) %>%
  summarise(position.lead=position.lead[1], W=sum(Z[!knockoff])-sum(Z[knockoff])) %>%
  ungroup() %>%
  arrange(desc(abs(W))) %>%
  mutate(Threshold = knockoff.threshold(W, fdr=0.1, offset=1)) %>%
  mutate(Selected = W >= Threshold) %>%
  mutate(causal=NA)

for(j in 1:nrow(Stats)) {
  bp <- Stats$position.lead[j]
  idx <- which((partition$bp.min<=bp)*(partition$bp.max>=bp)==1)
  Stats$causal[j] <- partition$causal[idx]
}

## Summarise discoveries
cat(sprintf("Preview of results:\n"))
Stats %>%
  filter(Selected) %>%
  summarise(N=n(), True=sum(causal), FDP=mean(!causal)) %>%
  print()

## Save lasso statistics
ofile <- sprintf("%s/a%s_s%s_%s_typed_%s_%s_lasso.txt",
                 odir, amplitude, seed, pop, density, res)
Lasso %>% write_delim(ofile, delim=" ")

## Save knockoff statistics
ofile <- sprintf("%s/a%s_s%s_%s_typed_%s_%s_stats.txt",
                 odir, amplitude, seed, pop, density, res)
Stats %>% write_delim(ofile, delim=" ")

cat(sprintf("Results written to:\n  %s/\n", odir))
