#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## Default parameters
phenotype <- "platelet"
resolution <- "res5"
q <- 0.1
r <- 2
nreps <- 2

## Input parameters
phenotype = as.character(args[1])
resolution = as.character(args[2])
q = as.numeric(args[3])
r = as.integer(args[4])
nreps <- as.integer(args[5])

set.seed(2021)

library(tidyverse)
library(gridExtra)
source("../paper/experiments/utils/functions_multienv.R")
source("../paper/experiments/utils/accumulation_test_functions.R")

## Load statistics
out.file <- "data/stats.txt"
if(FALSE) {
    ifiles.list <- Sys.glob("stats/lasso_*_*_related*.txt")
    data <- do.call("rbind", lapply(ifiles.list, function(ifile) {
        tmp <- strsplit(ifile, split="_")[[1]]
        resol <- strsplit(tmp[5], split=".txt")[[1]]
        df <- read_delim(sprintf("%s", ifile), delim=" ", col_types=cols()) %>%
            mutate(Phenotype = tmp[2], Population = tmp[3], Resolution = resol)
    }))
    ## Save combined statistics
    data %>% write_delim(out.file, delim=" ")
} else {
    ## Read combined statistics
    data <- read_delim(out.file, delim=" ", col_types=cols()) %>%
        filter(Phenotype == phenotype, Resolution == resolution)
}

## List of populations to process
pop.list <- c("asian", "black","britishindia","british","whitenonbritish")
stats <- data %>%
    filter(Population %in% pop.list) %>%
    group_by(Phenotype, Resolution, Population) %>%
    arrange(desc(abs(W))) %>%
    mutate(Threshold = knockoff.threshold(W, fdr=q, offset=1)) %>%
    mutate(Selected = W >= Threshold)
stats.everyone <- data %>%
    filter(Population == "everyone") %>%
    group_by(Phenotype, Resolution, Population) %>%
    arrange(desc(abs(W))) %>%
    mutate(Threshold = knockoff.threshold(W, fdr=q, offset=1)) %>%
    mutate(Selected = W > Threshold)

filter_partial_conjunction <- function(data, q, r=2, c=NULL, method="accumulation", randomness=3) {
    partial_conjunction_ <- function(stats, header, q, r=2, c=NULL, method="accumulation", randomness=3) {
        stats.tmp <- stats %>% arrange(CHR, Group) %>% mutate(CHR_Group = sprintf("%s_%s", CHR, Group))
        groups.list <- unique(stats.tmp$CHR_Group)
        n.groups <- length(groups.list)
        W.matrix <- matrix(0, n.groups, length(pop.list))
        for(k in 1:length(pop.list)) {
            pop <- pop.list[k]
            df <- filter(stats.tmp, Population==pop) %>% arrange(CHR, Group)
            groups.idx <- match(df$CHR_Group, groups.list)
            W.matrix[groups.idx,k] <- df$W
        }
        if(is.null(c)) {
            sel.pc <- partial_conjunction(W.matrix, r=r, q=q, method=method, randomness=randomness)
        } else {
            sel.pc <- partial_conjunction(W.matrix, r=r, q=q, c=c, randomness=randomness)
        }
        out <- stats.tmp %>% filter(CHR_Group %in% groups.list[sel.pc]) %>%
            group_by(CHR, Group) %>%
            summarise(idx=which.max(W), Population=Population[idx]) %>%
            ungroup()
        if(nrow(out)>1) {
            out <- as_tibble(cbind(header, out))
        } else {
            out <- header %>% mutate(CHR=NA, Group=NA, idx=NA, Population=as.character(NA)) %>% head(0)
        }
        return(out)
    }
    out <- data %>%
        group_by(Phenotype, Resolution) %>%
        group_map(~ partial_conjunction_(.x, .y, q, r=r, c=c, method=method, randomness=randomness), .keep = TRUE)
    out <- do.call("rbind", out)
    return(out)
}

stability_selection <- function(stats, q, r=2, c=NULL, method="accumulation", randomness=3, nreps=2) {
    discoveries <- do.call("rbind", lapply(1:nreps, function(rep) {
        cat(sprintf("Repetition %d of %d...\n", rep, nreps))
        df <- stats %>% filter_partial_conjunction(q, r=r, c=c, method=method, randomness=randomness)
        df %>% mutate(Rep = rep)
    }))
    results <- discoveries %>% group_by(Phenotype, Resolution, CHR, Group, Population)
}


## Collect results
res.num <- tibble()
res.disc <- tibble()

## Accumulation test (randomness 3)
df.tmp <- stats %>% stability_selection(q, r=r, method="accumulation", randomness=2, nreps=nreps)
num.tmp <- df.tmp %>% group_by(Phenotype, Resolution, Rep) %>% summarise(Discoveries=n(), nreps=nreps) %>%
    mutate(Method="Accumulation", Randomness=2, q=q, r=r)
disc.tmp <- df.tmp %>% group_by(Phenotype, Resolution, CHR, Group, Population) %>% summarise(Discovered = n(), nreps=nreps) %>%
    mutate(Method="Accumulation", Randomness=2, q=q, r=r)
res.num <- rbind(res.num, num.tmp)
res.disc <- rbind(res.disc, disc.tmp)

## Accumulation test (randomness 3)
df.tmp <- stats %>% stability_selection(q, r=r, method="accumulation", randomness=3, nreps=nreps)
num.tmp <- df.tmp %>% group_by(Phenotype, Resolution, Rep) %>% summarise(Discoveries=n(), nreps=nreps) %>%
    mutate(Method="Accumulation", Randomness=3, q=q, r=r)
disc.tmp <- df.tmp %>% group_by(Phenotype, Resolution, CHR, Group, Population) %>% summarise(Discovered = n(), nreps=nreps) %>%
    mutate(Method="Accumulation", Randomness=3, q=q, r=r)
res.num <- rbind(res.num, num.tmp)
res.disc <- rbind(res.disc, disc.tmp)

## SeqStep (randomness 2)
df.tmp <- stats %>% stability_selection(q, r=r, method="seqstep", c=0.5, randomness=2, nreps=nreps)
num.tmp <- df.tmp %>% group_by(Phenotype, Resolution, Rep) %>% summarise(Discoveries=n(), nreps=nreps) %>%
    mutate(Method="SeqStep", Randomness=2, q=q, r=r)
disc.tmp <- df.tmp %>% group_by(Phenotype, Resolution, CHR, Group, Population) %>% summarise(Discovered = n(), nreps=nreps) %>%
    mutate(Method="SeqStep", Randomness=2, q=q, r=r)
res.num <- rbind(res.num, num.tmp)
res.disc <- rbind(res.disc, disc.tmp)

## Intersection
disc.tmp <- stats %>% filter(Selected) %>%
    group_by(Phenotype, Resolution, CHR, Group) %>%
    summarise(N=n(), Population=paste(Population, collapse="_")) %>%
    ungroup() %>%
    filter(N>=r) %>%
    select(-N) %>%    
    mutate(Method="Intersection", Randomness=NA, q=q, r=r, nreps=1, Discovered=1)
res.disc <- rbind(res.disc, disc.tmp)

## Stack
disc.tmp <- stats.everyone %>% filter(Selected) %>%
    group_by(Phenotype, Resolution, CHR, Group) %>%
    mutate(Method="Stack", Randomness=NA, q=q, r=r, nreps=1, Discovered=1) %>%
    select(Phenotype, Resolution, CHR, Group, Population, Discovered, nreps, Method, Randomness, q, r)
res.disc <- rbind(res.disc, disc.tmp)

## Save results
out.file <- sprintf("results/discoveries_%s_%s_q%s_r%s.txt", phenotype, resolution, q, r)
res.disc %>% write_delim(out.file, delim=" ")
cat(sprintf("Results written to:\n  %s\n", out.file))
out.file <- sprintf("results/num_discoveries_%s_%s_q%s_r%s.txt", phenotype, resolution, q, r)
res.num %>% write_delim(out.file, delim=" ")
cat(sprintf("Results written to:\n  %s\n", out.file))
