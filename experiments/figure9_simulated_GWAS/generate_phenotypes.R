suppressMessages(library(tidyverse))
library(bigstatsr)

idir <- "/oak/stanford/groups/candes/matteo/invariant_knockoffs/data"
odir <- "/oak/stanford/groups/candes/matteo/invariant_knockoffs/data"

## Load sample file
ifile <- sprintf("%s/full_sample.txt", idir)
sample <- read_delim(ifile, delim=" ", col_types=cols())
populations <- unique(sample$population)
  
## Load legend file
ifile <- sprintf("%s/full_legend.txt", idir)
legend <- read_delim(ifile, delim=" ", col_types=cols())

## Load genotypes
ifile <- sprintf("%s/full_X.rds", idir)
fbm <- bigstatsr::big_attach(ifile)

#########################
## Generate phenotypes ##
#########################

source("utils.R")

## Load list of typed variants (highest density)
ifile <- sprintf("%s/typed_%s_legend.txt", idir, 0.1)
legend.typed <- read_delim(ifile, delim=" ", col_types=cols())
idx.typed <- which(legend$id %in% legend.typed$id)

## Number of causal variants (per population)
n.causal <- 10

## Choose causal variants
set.seed(2021)
legend.causal <- choose_causal(legend, idx.typed, populations, n.causal=n.causal)
legend.causal$sign <- (2*rbinom(nrow(legend.causal),1,0.5)-1)

## Estimate heritability
amplitude.list <- seq(10,100,by=10)
config <- expand.grid(amplitude=amplitude.list)
h2.table <- lapply(1:nrow(config), function(i) {
  seed <- config$seed[i]
  amplitude <- config$amplitude[i]
  tmp <- generate_phenotypes(fbm, legend.causal, amplitude=amplitude, seed=seed)
  out <- tibble(Amplitude=amplitude, h2 = tmp$h2)
  return(out)
})
h2.table <- do.call("rbind", h2.table)

## Store heritability table
ofile <- sprintf("%s/heritability.txt", odir)
h2.table %>% write_delim(ofile, delim=" ")
cat(sprintf("Heritability table written to:\n  %s\n", ofile))

## Generate phenotypes
seed.list <- 1:100
phenotypes <- tibble()
config <- expand.grid(amplitude=amplitude.list, seed=seed.list)
phenotypes <- lapply(1:nrow(config), function(i) {
  seed <- config$seed[i]
  amplitude <- config$amplitude[i]
  tmp <- generate_phenotypes(fbm, legend.causal, amplitude=amplitude, seed=seed)
  y <- tmp$y
  y <- round(as.matrix(y, length(y), 1),4)
  y.name <- sprintf("Y_a%d_s%d", amplitude, seed)
  colnames(y) <- y.name
  return(y)
})
phenotypes <- do.call("cbind", phenotypes)
phenotypes <- cbind(sample %>% filter(hap==0) %>% select(-hap), phenotypes) %>% as_tibble()

## Store phenotypes
ofile <- sprintf("%s/phenotypes.txt", odir)
phenotypes %>% write_delim(ofile, delim=" ")
cat(sprintf("Phenotypes written to:\n  %s\n", ofile))

## Store list of causal variants
ofile <- sprintf("%s/legend_causal.txt", odir)
legend.causal %>% write_delim(ofile, delim=" ")
cat(sprintf("Causal legend written to:\n  %s\n", ofile))

#######################
## Marginal analysis ##
#######################

if(FALSE) {

  y <- phenotypes[["Y_a100_s1"]]
  marginal.1 <- marginal_analysis(fbm, y, legend)

  marginal.1 %>%
    mutate(causal = id %in% legend.causal$id) %>%
    ggplot(aes(x=position/1e6, y=-log10(p.value), color=causal, shape=causal)) +
    geom_point(alpha=0.5, size=3) +
    geom_hline(yintercept=-log10(5e-8), linetype=2, color="blue") +
    scale_colour_manual(values = c("black", "red")) +
    theme_bw() +
    theme(text = element_text(size=20))

}

## marginal.1 <- marginal_analysis(data.typed, y, populations[1])
## marginal.2 <- marginal_analysis(data.typed, y, populations[2])

## pos.typed <- legend$position[idx.typed] / 1e6

## p1 <- marginal.1 %>%
##   mutate(position = pos.typed) %>%
##   ggplot(aes(x=position, y=-log10(p.value), color=causal, shape=causal)) +
##   geom_point(alpha=0.5, size=3) +
##   geom_hline(yintercept=-log10(5e-8), linetype=2, color="blue") +
##   scale_colour_manual(values = c("black", "red")) +
##   theme_bw() +
##   theme(text = element_text(size=20))

## p2 <- marginal.2 %>%
##   mutate(position = pos.typed) %>%
##   ggplot(aes(x=position, y=-log10(p.value), color=causal, shape=causal)) +
##   geom_point(alpha=0.5, size=3) +
##   geom_hline(yintercept=-log10(5e-8), linetype=2, color="blue") +
##   scale_colour_manual(values = c("black", "red")) +
##   theme_bw() +
##   theme(text = element_text(size=20))

## gridExtra::grid.arrange(p1, p2)
