## Experiments in Figure 6
## multi-environment knockoff filter: partial consistency testing
## FDR and Power vs amplitude
## Figure 6 runs the following parallelly with 
##              model in (1 4)
##              amplitude in (5 6.25 7.5 8.75 10 12.5 15 17.5 20 25 30)
##              main_seed in 1:100

main_seed = 1
amplitude = 20
model = 4
##################################################
library(tidyverse)
library(glmnet)
source("MEKF/functions_multienv.R")
source("MEKF/accumulation_test_functions.R")
source("MEKF/importance_stats.R")

plot.design <- function(beta, title, filename) {
    df <- as_tibble(beta) %>% mutate(Variable=1:nrow(beta)) %>%
        gather(V1, V2, V3, V4, V5, key="Environment", value="Effect") %>%
        mutate(Environment=str_replace(Environment, "V", ""))
    pp <- df %>%
        mutate(Effect=as.factor(Effect)) %>%
        ggplot(aes(x=Environment, y=Variable, fill=Effect, color=Effect)) +
        geom_tile() +
        scale_fill_manual(values = c("gray90", "black")) +
        scale_color_manual(values = c("gray90", "black")) +
        labs(color="Relative\nsignal\nstrength", fill="Relative\nsignal\nstrength") +
        theme_bw() +
        theme(legend.position = "right",
            plot.title = element_text(size = 10), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
            strip.text.x = element_text(size = 10), strip.text.y = element_text(size = 10), legend.title=element_text(size = 10)) +
        ggtitle(title)
    pp %>% ggsave(file=filename, width=3, height=2.25, units="in")
}

## Data parameters
q = 0.1
r = 2

if(model==4) {
    num_env = 5
    rhos = rep(0.2,num_env)
    n = 600
    p = 200
    beta = matrix(rep(0,num_env*p), nrow = p)
    k = 100
    beta[1:k,1:(num_env-1)] <- 1
    beta[seq(k+1,p),num_env] <- 1
    #beta[k+(1:k),num_env] <- 0.5
    #k0 = 40
    #for (i in 1:num_env){
    #    beta[2*k+(1:k0)+k0*(i-1), i] = 0.5
    #}
    ## Plot design
    if(FALSE) {
        plot.design(beta, "(b) Setting 2", sprintf("figures/beta_pinv_%s.pdf", model))
    }
} else {
    num_env = 5
    rhos = rep(0.2,num_env)
    n = 600
    p = 200
    beta = matrix(rep(0,num_env*p), nrow = p)
    k = 50
    k0 = 5
    beta[1:k,1:(num_env)] <- 1
    for (i in 1:num_env){
        beta[k+(1:k0)+k0*(i-1), i] = 1
    }
    ## Plot design
    if(FALSE) {
        plot.design(beta, "(a) Setting 1", sprintf("figures/beta_pinv_%s.pdf", model))
    }
}

beta = beta* amplitude / sqrt(n)

if(FALSE) {
    heatmap(beta,Rowv = NA,Colv = NA, scale = "none")
}

set.seed(2021)
beta.rows <- sample(1:nrow(beta))
beta <- beta[beta.rows,]
nonnulls = which(apply(beta, 1, combine_mag, r = r) > 0)

beta = beta* amplitude / sqrt(n)

out.file <- sprintf("results/m%s_a%d_s%d.txt", model, amplitude, main_seed)

## Define power and FDP
compute.power = function(selected, nonnulls) {mean(nonnulls %in% selected)}
compute.fdp = function(selected, nonnulls) {ifelse(length(selected)>0, 1-mean(selected %in% nonnulls), 0)}

run_experiment = function(seed) {
    set.seed(seed)

    ##------------------------
    ## Generate a dataset
    ##------------------------
    Xs = list()
    ys = list()
    for (i in 1:num_env){
        Sigma = toeplitz(rhos[i]^(0:(p-1)))
        X = matrix(rnorm(n*p),n) %*% chol(Sigma)
        Xs[[i]] = X
        invlogit = function(x) exp(x) / (1+exp(x))
        y.sample = function(x) rbinom(n, prob=invlogit(x %*% beta[,i]), size=1)
        y = y.sample(X)
        ys[[i]] = y
    }

    ##------------------------
    ## Create knockoffs
    ##------------------------
    X_ks = list()
    Sigma = toeplitz(rhos[1]^(0:(p-1)))
    diag_s = create.solve_sdp(Sigma)
    for (i in 1:num_env){
        X_ks[[i]] = create.gaussian(X, rep(0,p), Sigma, diag_s=diag_s)
    }

    ##------------------------
    ## Pool heuristic
    ##------------------------
    X_stack = NULL
    X_k_stack = NULL
    y_stack = NULL
    for (i in 1:num_env){
        X_stack = rbind(X_stack, Xs[[i]])
        X_k_stack = rbind(X_k_stack, X_ks[[i]])
        y_stack = c(y_stack,ys[[i]])
    }
    W_mag_all = stat.glmnet_coefdiff(X_stack, X_k_stack, y_stack, family="binomial")
    thres0 = knockoff.threshold(W_mag_all)
    selected = which(W_mag_all >= thres0)
    res.1 <- tibble(method="stack", power=compute.power(selected, nonnulls), fdp = compute.fdp(selected, nonnulls))

    # Compute stats environment by environment
    W_matrix = NULL
    for (i in 1:num_env){
        W = stat.glmnet_coefdiff(Xs[[i]], X_ks[[i]], ys[[i]], family="binomial")
        W_matrix = cbind(W_matrix,W)
    }

    ##------------------------
    ## Intersection heuristic
    ##------------------------
    selected01 = NULL
    for (i in 1:num_env){
        W_i = W_matrix[,i]
        thres_i = knockoff.threshold(W_i, fdr=q)
        selected01 = cbind(selected01, W_i >= thres_i)
    }
    selected = which(as.logical(apply(selected01, 1, combine_mag, r = r)))
    res.2 <- tibble(method="intersection", power=compute.power(selected, nonnulls), fdp = compute.fdp(selected, nonnulls))

    ##------------------------
    ## Partial consistency with meta-analysis stats
    ##------------------------
    selected = partial_conjunction(W_matrix, r = r, method = "accumulation")
    res.3 <- tibble(method="PI-meta-acc", power=compute.power(selected, nonnulls), fdp = compute.fdp(selected, nonnulls))
    selected = partial_conjunction(W_matrix, r = r, method = "accumulation", randomness = 3)
    res.4 <- tibble(method="PI-meta-acc-coin", power=compute.power(selected, nonnulls), fdp = compute.fdp(selected, nonnulls))
    selected = partial_conjunction(W_matrix, r = r, method = "seqstep", c=0.5)
    res.5 <- tibble(method="PI-meta-seqstep", power=compute.power(selected, nonnulls), fdp = compute.fdp(selected, nonnulls))

    ##------------------------
    ## Partial consistency with cross-prior stats
    ##------------------------
    W_cp_matrix = compute_stats_with_prior(Xs, X_ks, ys, verbose=TRUE, family="binomial")
    selected = partial_conjunction(W_cp_matrix, r = r, method = "accumulation", randomness = 3)
    res.6 <- tibble(method="PI-cp-acc", power=compute.power(selected, nonnulls), fdp = compute.fdp(selected, nonnulls))
    selected = partial_conjunction(W_cp_matrix, r = r, method = "accumulation")
    res.7 <- tibble(method="PI-cp-acc-coin", power=compute.power(selected, nonnulls), fdp = compute.fdp(selected, nonnulls))
    selected = partial_conjunction(W_cp_matrix, r = r, method = "seqstep", c=0.5)
    res.8 <- tibble(method="PI-cp-seqstep", power=compute.power(selected, nonnulls), fdp = compute.fdp(selected, nonnulls))

    ## Return results
    res = tibble(n=n, p=p, env=num_env, amplitude=amplitude, model=model, seed=seed) %>%
        cbind(rbind(res.1,res.2,res.3,res.4,res.5,res.6,res.7,res.8))

    return(res)
}

results = tibble()
for(num.exp in 1:1) {
    seed = main_seed*1000+num.exp
    res = run_experiment(seed)
    results = rbind(results, res)
    ## Save results
    results %>% write_delim(out.file, delim=" ")
    cat(sprintf("Results saved on %s\n", out.file))
}
