suppressMessages(library(tidyverse))

## Load legend
idir <- "/oak/stanford/groups/candes/matteo/invariant_knockoffs/data_1000gp/1000GP_Phase3"
ifile <- "1000GP_Phase3_chr22.legend.gz"
legend <- read_delim(sprintf("%s/%s", idir, ifile), delim=" ", col_types=cols())

## Load genetic map
ifile <- "genetic_map_chr22_combined_b37.txt"
map <- read_delim(sprintf("%s/%s", idir, ifile), delim=" ", col_types=cols())

## Interpolate missing genetic positions
interpolate.map <- approxfun(map$position, y=map$`Genetic_Map(cM)`)
legend$map <- interpolate.map(legend$position)
legend <- legend %>% select(id, position, map, everything())

## Filter variants
legend.filtered <- legend %>%
  filter(ALL > 0.001) %>% # Remove very rare variants
  filter(str_length(a0)==1, str_length(a1)==1) %>% # Keep only biallelic
  filter(TYPE=="Biallelic_SNP") %>%
  filter(!is.na(map)) ## Remove variants without map info

## Output list of variants to keep
odir <- "/oak/stanford/groups/candes/matteo/invariant_knockoffs/data_1000gp/1000GP_Phase3_clean"
ofile <- "legend_chr22.txt"
legend.filtered %>% write_delim(sprintf("%s/%s", odir, ofile), delim=" ")

idx.rows <- tibble(which(legend$id %in% legend.filtered$id))
ofile <- "rows_chr22.txt"
idx.rows %>% write_delim(sprintf("%s/%s", odir, ofile), delim=" ", col_names=FALSE)

## Extract appropriate rows from hap file
## sed -nf <(sed 's/$/p/' 1000GP_Phase3_clean/rows_chr22.txt) 1000GP_Phase3/1000GP_Phase3_chr22.hap | pv > 1000GP_Phase3_clean/1000GP_Phase3_chr22.hap

## Load sample file
ifile <- sprintf("%s/1000GP_Phase3.sample", idir, ifile)
sample <- read_delim(ifile, delim=" ", col_types=cols())

## Load haps file
idir <- "/oak/stanford/groups/candes/matteo/invariant_knockoffs/data_1000gp/1000GP_Phase3_clean"
ifile <- sprintf("%s/1000GP_Phase3_chr22.hap", idir, ifile)
df <- data.table::fread(ifile, sep=" ", showProgress=TRUE)
H <- as.matrix(df)
storage.mode(H) <- "integer"

## Convert haplotypes to FBM
ofile <- sprintf("%s/1000GP_Phase3_chr22", idir, ifile)
fbm <- bigstatsr::as_FBM(H, type = "integer", backingfile=ofile)
##bigstatsr::big_write(fbm, ofile, every_nrow = 100, progress = interactive())
fbm$save()
