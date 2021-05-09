## Experiments in Figure A11
## Transfer learning
## FDR and Power vs overlap
## Figure A11 runs the following with overlap = 0,20,40,60,80,100 parallelly

overlap = 100
##################################################
library(adaptiveKnockoff)
library(knockoff)
library(gam)
source("MEKF/importance_stats.R")
source("MEKF/multienv.R")

# Problem parameters
num_env = 3 ##the first env is the target env
ns = rep(800, num_env)
p = 500
q = 0.1
rhos = rep(0.5, num_env)
amplitude = 3.5
N = 500
k = 60
k_overlap = floor(k*overlap/100)

# Make the betas
beta = matrix(rep(0,num_env*p), nrow = p)
beta[1:k,1] = 1
beta[1:k + k - k_overlap,2:num_env] = 1
beta = beta* amplitude / sqrt(ns[1])
set.seed(1) ## make sure our model is not changing
permu = sample(1:p)
beta = beta[permu,]
nonnulls = which(beta[,1] != 0)

#heatmap(beta,Rowv = NA,Colv = NA, scale = "none")

fdps = NULL
powers = NULL
fdp = function(selected, nonnulls) length(setdiff(selected, nonnulls)) / max(1, length(selected))
power = function(selected, nonnulls) length(intersect(selected, nonnulls)) / length(nonnulls)                                                   

for (iter in 1:N){
    # generate data
    Xs = list()
    ys = list()
    for (i in 1:num_env){
        Sigma = toeplitz(rhos[i]^(0:(p-1)))
        X = matrix(rnorm(ns[i]*p),ns[i]) %*% chol(Sigma)
        Xs[[i]] = X
        y = X %*% beta[,i] + rnorm(ns[i])
        ys[[i]] = y
    }
    
    nonnulls = which(beta[,1]>0)
    
    ## create knockoffs
    X_ks = list()
    for (i in 1:num_env){
        X_ks[[i]] = create.gaussian(Xs[[i]], rep(0,p), toeplitz(rhos[i]^(0:(p-1))))
    }
    
    ## pool data from all environments
    X_stack = NULL
    X_k_stack = NULL
    y_stack = NULL
    for (i in 1:num_env){
        X_stack = rbind(X_stack, Xs[[i]])
        X_k_stack = rbind(X_k_stack, X_ks[[i]])
        y_stack = c(y_stack,ys[[i]])
    }
    W_stack = stat.lasso_coefdiff(X_stack, X_k_stack, y_stack)
    thres_stack = knockoff.threshold(W_stack, fdr = q)
    selected_stack = which(W_stack >= thres_stack)
    fdp_stack = fdp(selected_stack, nonnulls)
    power_stack = power(selected_stack,nonnulls)
    
    ## Vanilla knockoffs
    W_vanilla = stat.lasso_coefdiff(Xs[[1]], X_ks[[1]], ys[[1]])
    thres_vanilla = knockoff.threshold(W_vanilla, fdr = q)
    selected_vanilla = which(W_vanilla >= thres_vanilla)
    fdp_vanilla = fdp(selected_vanilla, nonnulls)
    power_vanilla = power(selected_vanilla,nonnulls)
    
    ## Our method -- transfer learning
    W_trans = compute_transfer_stats_with_prior(Xs, X_ks, ys)
    thres_trans = knockoff.threshold(W_trans, fdr = q)
    selected_trans = which(W_trans >= thres_trans)
    fdp_trans = fdp(selected_trans, nonnulls)
    power_trans = power(selected_trans,nonnulls)
    
    ## Our method -- MEKF
    W_i = compute_stats_with_prior(Xs, X_ks, ys)
    selected_i = invariant_model(W_i, q = q)
    fdp_i = fdp(selected_i, nonnulls)
    power_i = power(selected_i,nonnulls)
    
    ## Knockoffs with side information -- Zhimei Ren & Emmanuel Candes
    fdp_side = NA
    power_side = NA
    tryCatch({
        W_matrix_other = NULL
        for (i in 2:num_env){
            W = stat.lasso_coefdiff(Xs[[i]], X_ks[[i]], ys[[i]])
            W_matrix_other = cbind(W_matrix_other,W)
        }
        adap = filter_EM(W_vanilla, W_matrix_other, alpha = q)
        selected_side = adap$rejs[[1]]
        fdp_side = fdp(selected_side, nonnulls)
        power_side = power(selected_side, nonnulls)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
    fdps = cbind(fdps, c(fdp_stack, fdp_vanilla, fdp_trans, fdp_i, fdp_side))
    powers = cbind(powers, c(power_stack, power_vanilla, power_trans, power_i, power_side))
    #print(fdps)
    #print(powers)
}

fdps
powers

