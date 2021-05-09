library(tidyverse)
library(pbapply)
source("utils/functions_multienv.R")
source("utils/accumulation_test_functions.R")

## Load phenotype information
idir <- "/oak/stanford/groups/candes/matteo/invariant_knockoffs/data"
ifile <- sprintf("%s/heritability.txt", idir)
heritability <- read_delim(ifile, delim=" ", col_types=cols())

## Load partition information
ifile <- sprintf("%s/legend_partitions.txt", idir)
partitions <- read_delim(ifile, delim=" ", col_types=cols())

## Compute median group width in each partition
res.list <- 7:1
width.list <- sapply(res.list, function(res) {
  df <- partitions
  df$group <- partitions[[paste("res_",res, sep="")]]
  df <- df %>%
    group_by(group) %>%
    summarise(BP.min=min(position), BP.max=max(position), Width=(BP.max-BP.min)/1e3)
  return(round(median(df$Width)))
})
size.list <- sapply(res.list, function(res) {
  df <- partitions
  df$group <- partitions[[paste("res_",res, sep="")]]
  df <- df %>%
    group_by(group) %>%
    summarise(N=n())
  return(median(df$N))
})
df.partitions <- tibble(Resolution=res.list, Width=width.list, Size=size.list)

## Compute number of causal groups in each partition
idir <- "/oak/stanford/groups/candes/matteo/invariant_knockoffs/data"
ifile <- sprintf("%s/legend_causal.txt", idir)
causal <- read_delim(ifile, delim=" ", col_types=cols())
n.causal.list <- sapply(res.list, function(res) {
  df <- partitions
  df$group <- partitions[[paste("res_",res, sep="")]]
  df <- df %>% select(id, position, group) %>%
    left_join(causal %>% select(id, position) %>% mutate(causal=TRUE)) %>%
    mutate(causal = ifelse(is.na(causal), FALSE, causal)) %>%
    group_by(group) %>%
    summarise(causal=any(causal))
  n.causal <- sum(df$causal)
  return(n.causal)
})
df.n.causal <- tibble(Resolution=res.list, N.causal=n.causal.list)

#########################
## Load all statistics ##
#########################

idir <- "/oak/stanford/groups/candes/matteo/invariant_knockoffs/stats"
odir <- "/oak/stanford/groups/candes/matteo/invariant_knockoffs/results"

if(FALSE) {  

  library(parallel)
  cl <- makeCluster(10)
  
  ifiles.list <- list.files(idir, pattern = "stats.txt$") 
  clusterExport(cl, c("idir", "ifiles.list", "tibble", "parse_number", "read_delim", "cols", "col_integer", "col_double", "col_logical"))

  data <- do.call("rbind", pblapply(1:length(ifiles.list), cl=cl, function(i) {
    ifile <- ifiles.list[i]
    x <- strsplit(ifile, "_")[[1]]
    header <- tibble(Amplitude=parse_number(x[1]),
                     Seed=parse_number(x[2]),
                     Density=parse_number(x[5]),
                     Resolution = parse_number(x[7]),
                     Population=x[3]
                     )
    df <- read_delim(sprintf("%s/%s", idir, ifile), delim=" ",
                     col_types=cols(
                       group = col_integer(),
                       position.lead = col_integer(),
                       W = col_double(),
                       Threshold = col_double(),
                       Selected = col_logical(),
                       causal = col_logical()
                     ))
    df <- cbind(header, df)
    return(df)
  }))
  data <- as_tibble(data)
  
  data <- data %>% left_join(df.n.causal) %>% mutate(Resolution = factor(Resolution, res.list, width.list))

  data <- data %>% left_join(heritability) %>% select(Amplitude, h2, everything())  
  ## Save results
  ofile <- sprintf("%s/stats.txt", odir)
  data %>% write_delim(ofile, delim=" ")

} else {
  ofile <- sprintf("%s/stats.txt", odir)
  data <- read_delim(ofile, delim=" ")  
}

header <- data %>%
  group_by(Amplitude, h2, Seed, Density, Resolution) %>%
  summarise()

## Apply knockoff filter
q <- 0.2

data <- data %>%
  group_by(Amplitude, h2, Seed, Density, Resolution, N.causal, Population) %>%
  arrange(desc(abs(W))) %>%
  mutate(Threshold = knockoff.threshold(W, fdr=q, offset=1)) %>%
  mutate(Selected = W >= Threshold)


########################
## Invariant analysis ##
########################

filter_partial_conjunction <- function(data, q, r=2, c=NULL, method="accumulation", randomness=2) {
  partial_conjunction_ <- function(stats, header, q, r=2, c=NULL, method="accumulation", randomness=2) {
    n.groups <- max(stats$group)
    W.matrix <- matrix(0, n.groups, length(pop.list))
    for(k in 1:length(pop.list)) {
      pop <- pop.list[k]
      df <- filter(stats, Population==pop) %>% arrange(group)
      groups <- df$group
      W.matrix[groups,k] <- df$W
    }
    if(is.null(c)) {
      sel.pc <- partial_conjunction(W.matrix, r=r, q=q, method=method, randomness=randomness)
    } else {
      sel.pc <- partial_conjunction(W.matrix, r=r, q=q, c=c, randomness=randomness)
    }
    out <- stats %>% filter(group %in% sel.pc) %>%
      group_by(group, causal) %>%
      summarise(idx=which.max(W), Population=Population[idx]) %>%
      ungroup()
    if(nrow(out)>1) {
      out <- as_tibble(cbind(header, out))
    } else {
      out <- header %>% mutate(group=NA, causal=NA, idx=NA, Population=NA) %>% head(0)
    }
    return(out)
  }
  out <- data %>%
    group_by(Amplitude, h2, Seed, Density, Resolution, N.causal) %>%
    group_map(~ partial_conjunction_(.x, .y, q, r=r, c=c, method=method), .keep = TRUE)
  out <- do.call("rbind", out)
  return(out)
}

r <- 3
pop.list <- c("EUR", "EAS", "AMR", "SAS", "AFR")

## Partial conjunction
#disc.pc <- data %>% filter_partial_conjunction(q, r=r, method="accumulation")
disc.pc <- data %>%
#  filter(Density==0.01, Resolution%in%c(1,35,463), Amplitude %in% c(10,30,50,80,100)) %>%
  filter_partial_conjunction(q, r=r, method="accumulation", randomness=2)
#  filter_partial_conjunction(q, r=r, method="seqstep", c=0.75)
res.pc <- disc.pc %>%
  group_by(Amplitude, h2, Seed, Density, Resolution, N.causal) %>%
  summarise(N=n(), True=sum(causal), FDP=mean(!causal)) %>%
  ungroup()
res.pc <- header %>% left_join(res.pc) %>%
  replace_na(list(N = 0, True = 0, FDP=0)) %>%
  mutate(Method = "Partial conjunction")

## Merged analysis
res.merged <- data %>% filter(Population=="merged", Selected) %>%
  group_by(Amplitude, h2, Seed, Density, Resolution, N.causal) %>%
  summarise(N=n(), True=sum(causal), FDP=mean(!causal)) %>%
  ungroup()
res.merged <- header %>% left_join(res.merged) %>%
  replace_na(list(N = 0, True = 0, FDP=0)) %>%
  mutate(Method = "Stack")
  
## Selections in at least r different environments
res.int <- data %>%
  filter(Selected) %>%
  group_by(Amplitude, h2, Seed, Density, Resolution, N.causal, group, causal) %>%
  summarise(N=n(), Population=paste(Population, collapse="_")) %>%
  ungroup() %>%
  filter(N>=r) %>%
  group_by(Amplitude, h2, Seed, Density, Resolution, N.causal) %>%
  summarise(N=n(), True=sum(causal), FDP=mean(!causal)) %>%
  ungroup()
res.int <- header %>% left_join(res.int) %>%
  replace_na(list(N = 0, True = 0, FDP=0)) %>%
  mutate(Method = "Intersection")

## Save results
results.joint <- rbind(res.pc, res.merged, res.int)
ofile <- sprintf("%s/results_joint.txt", odir)
results.joint %>% write_delim(ofile, delim=" ")

################
## Make plots ##
################

rbind(res.pc, res.merged, res.int) %>%
  group_by(Amplitude, h2, Density, Resolution, Method) %>%
  summarise(N=mean(N), True=mean(True), FDP=mean(FDP)) %>%
  ggplot(aes(x=h2, y=FDP, color=Method)) +
  geom_point() +
  geom_line() +
  facet_grid(Density~Resolution) +
  geom_hline(yintercept=q) +
  scale_x_sqrt() +
  ylim(0,1) +
  theme_bw()

rbind(res.pc, res.merged, res.int) %>%
  mutate(N.causal=ifelse(is.na(N.causal),1,N.causal), Power=True/N.causal) %>%  
  group_by(Amplitude, h2, Density, Resolution, Method) %>%
  summarise(N=mean(N), True=mean(True), Power=mean(Power), FDP=mean(FDP)) %>%
  ggplot(aes(x=h2, y=Power, color=Method)) +
  geom_point() +
  geom_line() +
  scale_x_sqrt() +
  facet_grid(Density~Resolution) +
  theme_bw()


##############
## Misc     ##
##############

## Selections in each environment separately
res.sep <- data %>%
  filter(Selected) %>%
  group_by(Amplitude, h2, Seed, Density, Resolution, Population, N.causal) %>%
  summarise(N=n(), True=sum(causal), FDP=mean(!causal)) %>%
  mutate(N.causal=ifelse(is.na(N.causal),1,N.causal), Power=True/N.causal)
  

## Save results
ofile <- sprintf("%s/results_separate.txt", odir)
res.sep %>% write_delim(ofile, delim=" ")

res.sep %>%
  group_by(Amplitude, h2, Density, Resolution, Population) %>%
  summarise(N=mean(N), True=mean(True), FDP=mean(FDP)) %>%
  ggplot(aes(x=h2, y=FDP, color=Population)) +
  geom_point() +
  geom_line() +
  facet_grid(Density~Resolution) +
  geom_hline(yintercept=q) +
  ylim(0,1) +
  theme_bw()

res.sep %>%
  group_by(Amplitude, h2, Density, Resolution, Population) %>%
  summarise(N=mean(N), True=mean(True), FDP=mean(FDP)) %>%
  ggplot(aes(x=h2, y=True, color=Population)) +
  geom_point() +
  geom_line() +
  facet_grid(Density~Resolution) +
  theme_bw()
