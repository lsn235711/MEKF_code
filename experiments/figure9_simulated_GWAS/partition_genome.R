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

## Values of resolution (measured in cM)
resolution.list <- c(0, 0.01, 0.05, 0.1, 0.2, 0.5, 1, 2)

legend <- legend %>%
  mutate(MAF=pmin(ALL,1-ALL)) %>%
  select(id, position, map, MAF)

## Divide variants into blocks by K-means clustering
n.clusters <- ceiling(nrow(legend)/10000)
clustering <- kmeans(legend$map, n.clusters, nstart = 20)
clusters <- factor(clustering$cluster, levels=unique(clustering$cluster), labels=1:n.clusters) %>%
  as.integer()
cluster.sizes <- sapply(1:max(clusters), function(c) sum(clusters==c))

## Hierarchical clustering within each block
hc.list <- lapply(1:n.clusters, function(c) {
  c.idx <- which(clusters==c)
  map <- legend$map[c.idx]
  D <- dist(map)
  hc <- fastcluster::hclust(D, method="complete")
  return(hc)
})

cut.trees <- function(hc.list, h) {
  ## Cut each tree at a certain height
  hc.clusters.list <- lapply(1:n.clusters, function(c) {
    hc.clusters <- cutree(hc.list[[c]], h=h)
    return(hc.clusters)
  })
  ## Combine list of clusters
  partial.sum <- 0
  hc.clusters <- c()
  for(c in 1:n.clusters) {
    hc.clusters <- c(hc.clusters, hc.clusters.list[[c]] + partial.sum)
    partial.sum <- partial.sum + max(hc.clusters.list[[c]])
  }
  return(hc.clusters)
}

## Cut dendrogram to define groups
Partitions <- legend %>% select(id, position, map)
for(resolution.idx in 1:length(resolution.list)) {
  resolution <- resolution.list[resolution.idx]
  col.name <- sprintf("res_%d", length(resolution.list)-resolution.idx+1)
  hc.clusters <- cut.trees(hc.list, resolution)    
  Partitions[[col.name]] <- hc.clusters
}

## Compute average group sizes
group.sizes <- apply(Partitions[,-c(1,2,3)], 2, function(x) mean(table(x)))
cat(sprintf("Mean group sizes: \n"))
print(group.sizes)

## Save results
ofile <- sprintf("%s/legend_partitions.txt", odir)
Partitions %>% write_delim(ofile, delim=" ")
cat(sprintf("Partitions written to: %s\n", ofile))
