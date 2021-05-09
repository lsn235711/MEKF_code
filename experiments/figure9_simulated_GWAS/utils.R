library(tidyverse)

generate_phenotypes <- function(fbm, legend.causal, amplitude=20, seed=1) {
  set.seed(seed)
  n <- nrow(fbm)
  beta <- rep(amplitude, nrow(legend.causal)) * legend.causal$sign / sqrt(n)
  idx.causal <- which(legend$id %in% legend.causal$id)
  X <- fbm[,idx.causal]
  storage.mode(X) <- "double"
  X <- scale(X)
  if(sum(is.na(X))>0) {
    cat("Warning: found NA")
    X[is.na(X)] <- 0
  }
  y.mean <- X %*% beta
  ## Estimate heritability and choose interesting signal amplitudes
  y.mean.var <- as.numeric(var(y.mean))
  amplitude <- mean(abs(beta))
  h2 <- amplitude^2*y.mean.var / (amplitude^2*y.mean.var + 1)
  ## Generate noise
  y <- y.mean + rnorm(n)
  ## Return phenotype and heritability
  out <- c()
  out$y <- y
  out$h2 <- h2
  return(out)
}

compute_typing_indices <- function(legend, density.list, seed=2021) {
  set.seed(seed)

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

  prototypes.list <- lapply(density.list, function(density) {
    prototypes.list <- lapply(1:n.clusters, function(c) {
      ## Cut each tree to obtain the desired number of groups
      k <- ceiling(cluster.sizes[c]*density)
      hc.clusters <- cutree(hc.list[[c]], k=k)
      ## Pick one representative per group: maximize MAF
      n.groups <- max(hc.clusters)
      prototypes.tmp <- sapply(1:n.groups, function(i) {
        idx.group <- which(hc.clusters==i)
        j <- idx.group[which.max(legend$MAF[idx.group])]
      })
      return(prototypes.tmp)
    })
    ## Combine list of prototypes
    partial.sum <- 0
    prototypes <- c()
    for(c in 1:n.clusters) {
      prototypes <- c(prototypes, prototypes.list[[c]] + partial.sum)
      partial.sum <- partial.sum + cluster.sizes[c]
    }
    return(prototypes)
  })

  return(prototypes.list)
}

choose_causal <- function(legend, idx.typed, populations, n.causal=10) {
  p <- nrow(legend)
  idx.missing <- setdiff(1:p, idx.typed)
  legend.miss <- legend[idx.missing,]
  df <- as.matrix(legend.miss[,populations])
  legend.miss$pop <- populations[apply(df,1,which.max)]
  legend.miss$maf <- apply(df,1,max)
  legend.miss$maf <- pmin(legend.miss$maf, 1-legend.miss$maf)
  legend.miss <- legend.miss %>% select(id, position, map, pop, maf, everything())
  maf.min <- 0.1
  maf.max <- 0.5
  legend.causal <- legend.miss %>%
    filter(maf>=maf.min, maf<=maf.max) %>%
    group_by(pop) %>%
    sample_n(n.causal) %>%
    ungroup() %>%
    arrange(position)
  return(legend.causal)
}

extract_typed_data <- function(data.full, idx.typed) {
  p <- ncol(data.full$X)
  data.typed <- c()
  data.typed$Z <- data.full$Z[,idx.typed]
  data.typed$H <- data.full$H[,idx.typed]
  data.typed$X <- data.full$X[,idx.typed]
  data.typed$population <- data.full$population
  return(data.typed)
}

generate_haplotypes_onepop <- function(hmm, idx.row) {
  n <- length(idx.row)
  p <- nrow(hmm$alpha)
  K <- ncol(hmm$alpha)
  alpha <- hmm$alpha[1,]
  ## Sample the haplotypes
  ## Transition probabilities
  ej <- exp(-hmm$r)
  prob.jump <- c(1, 1-(ej + (1-ej)/K))
  ## Decide where to carry out random transitions
  z <- rep(1,n)
  pb <- txtProgressBar(min = 0, max = p, style = 3, file=stderr())
  idx.row.x <- round(idx.row[seq(2,length(idx.row),by=2)]/2)
  indices.h1 <- seq(1,length(idx.row),by=2)
  indices.h2 <- seq(2,length(idx.row),by=2)
  for(j in 1:p) {
    jump <- rbinom(n, 1, prob.jump[j])
    n.jumps <- sum(jump)
    z[jump==1] <- sample(1:K, n.jumps, replace=TRUE, prob=alpha)
    h <- rbinom(n, 1, hmm$theta[j,z])
    H[idx.row,j] <- h
    X[idx.row.x,j] <- h[indices.h1] + h[indices.h2]
    setTxtProgressBar(pb, j)
  }
  close(pb)
  return(NULL)
}

generate_genotypes <- function(hmm.list, legend, n, return.Z=FALSE) {
  K <- length(hmm.list)
  pop.idx <- rep(1:K, length.out=n)
  pop.n <- sapply(1:K, function(k) sum(pop.idx==k))
  out <- c()
  out$H <- lapply(1:K, function(k) generate_haplotypes_onepop(hmm.list[[k]], pop.n[k]))
  out$H <- do.call("rbind", out$H)
  out$X <- out$H[seq(1,2*n,by=2),] + out$H[seq(2,2*n,by=2),]
  out$legend <- legend
  population <- do.call("c", (lapply(1:K, function(k) rep(hmm.list[[k]]$population, 2*pop.n[k]))))
  out$sample <- tibble(id=rep(1:n, each=2), hap=rep(c(0,1), length.out=nrow(out$H)), population=population)
  return(out)
}

initialize_hmm_single <- function(fbm.ref, legend, legend.full) {
  K <- nrow(fbm.ref)
  p.full <- ncol(fbm.ref)
  p <- nrow(legend)
  idx.typed <- which(legend.full$id %in% legend$id)
  r.scale <- 1
  mutation.rate <- 1e-3
  out <- c()
  out$alpha <- array(1/K, c(p, K))
  out$r <- c(0,r.scale * diff(legend$map))
  out$theta <- t(fbm.ref[,idx.typed])
  out$theta <- pmax(pmin(out$theta,1-mutation.rate),mutation.rate)
  return(out)
}

initialize_hmm <- function(ref1000gp) {
  K <- length(ref1000gp)
  p <- ncol(ref1000gp[[1]]$H)
  M <- nrow(ref1000gp[[1]]$H)

  r.scale <- 1
  mutation.rate <- 1e-3

  hmm <- lapply(1:K, function(k) {
    out <- c()
    out$alpha <- array(1/M, c(p, M))
    out$r <- r.scale * diff(ref1000gp[[k]]$variants$map)
    out$theta <- t(ref1000gp[[k]]$H)
    out$theta <- pmax(pmin(out$theta,1-mutation.rate),mutation.rate)
    out$population <- ref1000gp[[k]]$population
    return(out)
  })

  return(hmm)
}

subset_1000gp <- function(fbm, sample, legend, populations, M=2, p=NULL) {
  if(is.null(p)) p <- nrow(legend)

  lapply(populations, function(population) {
    ## Initialize output container
    dat.out <- c()

    ## Copy list of variants
    dat.out$variants <- legend[1:p,]

    ## Extract list of samples
    sample.tmp <- sample %>%
      filter(GROUP==population) %>%
      sample_n(M)
    dat.out$sample <- sample %>%
      inner_join(sample.tmp, by = colnames(sample))

    ## Store population information
    dat.out$population <- population

    ## Copy haplotype data
    idx.col <- which(sample$ID %in% dat.out$sample$ID)
    dat.out$H <- t(fbm[1:p,idx.col])

    return(dat.out)
  })
}

marginal_analysis <- function(fbm, y, legend, ind.row=NULL) {
  if(is.null(ind.row)) ind.row <- 1:length(y)
  marginal <- bigstatsr::big_univLinReg(fbm, y, ind.train=ind.row)
  out <- legend
  out$p.value <- predict(marginal, log10 = FALSE)
  out <- as_tibble(out) %>%
    select(id, position, map, p.value, everything())
  return(out)
}

transform_hmm <- function(hmm, indices) {
  K <- length(hmm$pInit)
  indices <- sort(indices)
  hmm_new <- c()
  pInit.tmp <- hmm$pInit
  if(indices[1] > 1) {
    for(j in 1:(indices[1]-1)) {
      pInit.tmp <- t(hmm$Q[j,,]) %*% pInit.tmp
    }
  }
  hmm_new$pInit <- as.numeric(pInit.tmp)

  hmm_new$Q <- array(NA,c(length(indices)-1,K,K))
  for(j in 1:(length(indices)-1)) {
    Q.prod <- diag(length(hmm$pInit))
    if(j==1){
      start = 1
    } else {
      start = indices[j-1]
    }
    for(k in start:(indices[j]-1)) {
      Q.prod <- Q.prod %*% hmm$Q[j,,]
    }
    hmm_new$Q[j,,] = Q.prod
  }
  ## Emission distribution
  hmm_new$pEmit <- hmm$pEmit[indices,]
  ## Population
  hmm_new$populations <- hmm$populations
  return(hmm_new)
}

visualize_haplotypes <- function(data, hmm) {
  K <- length(hmm$populations)
  M <- ncol(hmm$Q)/K

  ## Define the color ramp (returns a function object)
  ramp <- colorRamp(c("blue", "red"))

  ## Define the ramp hex triplets
  col.list <- rgb( ramp(seq(0, 1, length = K)), max = 255)

  superheat::superheat(t(hmm$pEmit), scale=F, heat.col.scheme="grey",
                       row.dendrogram = TRUE,
                       order.rows=1:ncol(hmm$pEmit),
                       legend = FALSE,
                       left.label.col=rep(col.list, each=M),
                       left.label.text.size = 3.5)

  package.trajectories <- function(Z, pop) {
    i <- 1
    df1 <- tibble(Sample=i, Variable = 1:ncol(Z), State=Z[i,], Population=pop[i])
    i <- 2
    df2 <- tibble(Sample=i, Variable = 1:ncol(Z), State=Z[i,], Population=pop[i])
    i <- 3
    df3 <- tibble(Sample=3, Variable = 1:ncol(Z), State=Z[i,], Population=pop[i])
    df <- rbind(df1, df2, df3)
  }

  df <- package.trajectories(data$Z, data$population) %>% mutate(Model="Pop. structure")
  p.trajectories <- df %>%
    ggplot(aes(x=Variable, y=State, color=Population, group=1)) +
    geom_point() + geom_line() +
    facet_grid(Model~Sample, labeller = labeller(Sample=label_both), switch="y") +
    scale_colour_manual(values = col.list) +
    ylim(1,M*K) +
    scale_x_continuous(breaks = c(0,250,500)) +
    ylab("Markov chain state") +
    theme_bw() +
    theme(legend.position="right", text = element_text(size=25))
  p.trajectories
}
