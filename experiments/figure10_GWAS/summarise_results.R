library(tidyverse)
library(gridExtra)
library(kableExtra)

if(FALSE) {
    ifile.list <- Sys.glob("results/discoveries_*.txt")
    results <- do.call("rbind", lapply(ifile.list, function(ifile) {
        df <- read_delim(sprintf("%s", ifile), delim=" ", col_types=cols())
    }))
    ## Load partitions
    ifile <- "partitions.txt"
    Partitions <- read_delim(ifile, delim=" ") %>% mutate(Resolution = sprintf("res%d", Resolution))
    Groups <- Partitions %>% group_by(Resolution, CHR, Group) %>% summarise(BP.min=min(BP), BP.max=max(BP)) %>% ungroup()
    ## Cross-references discoveries with groups
    results <- results %>% left_join(Groups, by = c("Resolution", "CHR", "Group"))
    res.file <- "discoveries.txt"
    results %>% write_delim(res.file, delim=" ")
    ifile.list <- Sys.glob("results/num_discoveries_*.txt")
    num.results <- do.call("rbind", lapply(ifile.list, function(ifile) {
        df <- read_delim(sprintf("%s", ifile), delim=" ", col_types=cols())
    }))
    res.file <- "num_discoveries.txt"
    num.results %>% write_delim(res.file, delim=" ")
} else {
    res.file <- "discoveries.txt"
    results <- read_delim(res.file, delim=" ", col_types=cols())
    res.file <- "num_discoveries.txt"
    num.results <- read_delim(res.file, delim=" ", col_types=cols(),  guess_max=10000) %>%
        mutate(Discoveries = ifelse(is.na(Discoveries), N, Discoveries)) %>%
        select(-N, -nreps)
}

df <- results %>%
    filter(Discovered > nreps / 2) %>%
    group_by(Phenotype, Resolution, Method, Randomness, q, r) %>%
    summarise(Discoveries = n())

## df %>%
##     filter(r==2, Method!="Stack", ! Randomness %in% c(3)) %>%
##     ggplot(aes(x=Method, y=Discoveries)) +
##     geom_col() +
##     facet_grid(Phenotype ~ Resolution, scales="free") +
##     theme_bw()

res.values <- c("res6", "res5", "res4", "res3", "res2", "res1", "res0") %>% rev()
res.labels <- c(425, 208, 81, 41, 20, 3, "single-SNP") %>% rev()
r.values <- 1:5

out.dir <- "../paper/manuscript/tables"
out.dir.fig <- "../paper/manuscript/figures_matteo"

######################################################
## Table 1: Discoveries for height and platelet vs resolution, r
######################################################

nominal_fdr <- 0.1

df.stack <- df %>%
    filter(Method=="Stack", q==nominal_fdr, r==5) %>%
    ungroup() %>%
    select(-Randomness, -q, -Method) %>%
    mutate(Resolution = factor(Resolution, res.values, res.labels), r=1) %>%
    arrange(Phenotype, Resolution, r) %>%
    complete(Phenotype, Resolution) %>%
    replace_na(list(Discoveries = 0, r=1))

df1 <- df %>%
    filter(Method=="Accumulation", Randomness==2, q==nominal_fdr) %>%
    ungroup() %>%
    select(-Randomness, -q, -Method) %>%
    mutate(Resolution = factor(Resolution, res.values, res.labels)) %>%
    rbind(df.stack) %>%
    mutate(r = factor(r, r.values, r.values)) %>%
    arrange(Phenotype, r, Resolution) %>%
    complete(Phenotype, r, Resolution, ) %>%
    replace_na(list(Discoveries = 0))

tb <- df1 %>%
    spread(r, Discoveries)

tbl <- tb %>%
    filter(Phenotype %in% c("platelet", "height")) %>%
    kbl("latex", align="c", booktabs=TRUE, escape=T,
                  col.names = c("Phenotype", "Resolution (kb)", 1:5)) %>%
    add_header_above(c(" " = 2, "Number of environments" = 5)) %>%
    collapse_rows(columns = 1, valign = "top", latex_hline="major", longtable_clean_cut=TRUE) %>%
    column_spec(column = 3:7, width = "0.7cm", latex_valign = "m")


out.file <- sprintf("%s/analysis_discoveries_small.tex", out.dir)
tbl %>% save_kable(out.file, keep_tex=TRUE, self_contained=FALSE)


######################################################
## Table 2: like Table 1, with all phenotypes
######################################################

tbl <- tb %>%
    filter(!Phenotype %in% c("platelet", "height")) %>%
    kbl("latex", align="c", booktabs=TRUE, escape=T,
                  col.names = c("Phenotype", "Resolution (kb)", 1:5)) %>%
    add_header_above(c(" " = 2, "Number of environments" = 5)) %>%
    collapse_rows(columns = 1, valign = "top", latex_hline="major", longtable_clean_cut=TRUE) %>%
    column_spec(column = 3:7, width = "0.7cm", latex_valign = "m")


out.file <- sprintf("%s/analysis_discoveries.tex", out.dir)
tbl %>% save_kable(out.file, keep_tex=TRUE, self_contained=FALSE)

######################################################
## Table 3: Discoveries for height and platelet vs resolution, r, other heuristics
######################################################

r.values <- 2:5
method.values <- c("Accumulation", "SeqStep", "Intersection")
method.labels <- c("Accumulation", "SeqStep", "Intersection")
method.labels.short <- c("Acc.", "SStep", "Int.")

df0 <- expand.grid(r=c(2:5), Method=method.values) %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    arrange(r, Method) %>%
    mutate(Index = sprintf("%s_%s", r, Method)) %>% select(-Method, -r)
index.values <- df0$Index
df0

df1 <- df %>%
    filter(!Randomness %in% c(3), Method != "Stack", q==nominal_fdr) %>%
    ungroup() %>%
    select(-Randomness, -q) %>%
    mutate(Resolution = factor(Resolution, res.values, res.labels)) %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    mutate(r = factor(r, r.values, r.values)) %>%
    arrange(Phenotype, r, Resolution, Method) %>%
    complete(Phenotype, r, Resolution, Method) %>%
    replace_na(list(Discoveries = 0)) %>%
    mutate(Index = sprintf("%s_%s", r, Method)) %>%
    select(-r, -Method)

tb <- df1 %>%
    spread(Index, Discoveries) %>%
    select(Phenotype, Resolution, df0$Index)

tbl <- tb %>%
    kbl("latex", align="c", booktabs=TRUE, escape=T,
                  col.names = c("Phenotype", "Resolution (kb)", rep(method.labels.short, 4))) %>%
    add_header_above(c(" " = 2, "2" = 3, "3" = 3, "4" = 3, "5" = 3)) %>%
    add_header_above(c(" " = 2, "Number of environments / Invariance testing method" = 12)) %>% 
    collapse_rows(columns = 1, valign = "top", latex_hline="major", longtable_clean_cut=TRUE) %>%
    column_spec(column = 3:7, width = "0.5cm", latex_valign = "m")

out.file <- sprintf("%s/analysis_discoveries_methods.tex", out.dir)
tbl %>% save_kable(out.file, keep_tex=TRUE, self_contained=FALSE)


######################
## Plot variability ##
######################

nominal_fdr <- 0.1
r_value <- 2
res_value <- "res6"

method.values <- c("Accumulation - 3", "Accumulation - 2", "SeqStep - 2", "Intersection - 0")
method.labels <- c("Accumulation (coin)", "Accumulation", "SeqStep", "Intersection")

df <- results %>%
    filter(Discovered > nreps / 2) %>%
    group_by(Phenotype, Resolution, Method, Randomness, q, r) %>%
    summarise(Discoveries = n())
num.int <- df %>%
    filter(Method=="Intersection", q==nominal_fdr, r==r_value) %>%
    ungroup() %>%
    arrange(Phenotype, Resolution, r) %>%
    complete(Phenotype, Resolution) %>%
    replace_na(list(Discoveries = 0, r=r_value, Randomness = 0)) %>%
    filter(Resolution==res_value) %>%
    mutate(Rep=0)

pp <- num.results %>%
    mutate(Phenotype = ifelse(Phenotype=="hypothyroidism", "hypothyr.", Phenotype)) %>%
    rbind(num.int) %>%
    filter(Resolution=="res6", Method %in% c("Accumulation", "SeqStep", "Intersection"), q==0.1, r==2) %>%
    mutate(Method = sprintf("%s - %s", Method, Randomness)) %>%
    complete(Phenotype, Method) %>%
    mutate(Method = factor(Method, method.values, method.labels)) %>%
    replace_na(list(Discoveries = 0)) %>%
    ggplot(aes(x=Method, y=Discoveries)) +
    geom_boxplot() +
    facet_wrap(.~Phenotype, scales="free", nrow=1) +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(size = 10), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
          axis.text.x = element_text(angle = 90),
          strip.text.x = element_text(size = 10), strip.text.y = element_text(size = 10), legend.title=element_text(size = 10))
   
pp %>% ggsave(file=sprintf("%s/analysis_numdisc.pdf", out.dir.fig), width=9, height=3.5, units="in")


##################
## Chicago plot ##
##################

source("utils_chicagoplot.R")

df.1 <- results %>%
    filter((Discovered > 50) | (Method=="Stack")) %>%
    mutate(r = ifelse(Method=="Stack", 1, r)) %>%
    mutate(Randomness = ifelse(Method=="Stack", 2, Randomness)) %>%
    filter(Phenotype=="platelet", Method %in% c("Accumulation", "Stack"), Randomness==2, q==0.1) %>%
    select(-Phenotype, -q) %>%
    mutate(Method="Knockoffs", FDP.local=0, Resolution=parse_number(Resolution))

window.chr <- 3
window.left <-55.5e6
window.right <-59e6
p.knockoffs.overlay <- plot_chicago(window.chr, window.left, window.right, df.1, upside=TRUE, plot.title="", overlay=TRUE)
p.knockoffs.overlay %>% ggsave(file="platelet_chicago.pdf", width=7, height=2.25, units="in")

p.knockoffs <- plot_chicago(window.chr, window.left, window.right, df.1, upside=TRUE, plot.title="", overlay=FALSE)
p.knockoffs %>% ggsave(file="platelet_chicago_stack.pdf", width=6, height=5.5, units="in")
