## Experiments in Figure A2
## An example of invalid multi-environment statistics
## FDR and Power vs amplitude
## Figure A2 runs the following with amp_no = 1:10 parallelly

amp_no = 8
##################################################
library(knockoff)
source("MEKF/functions_multienv.R")
source("MEKF/accumulation_test_functions.R")
source("MEKF/importance_stats.R")

# Problem parameters
num_env = 2 ##the first env is the target env
ns = rep(200,2)
p = 100
q = 0.1
rhos = rep(0.6,num_env)
amplitude = amp_no

N = 500
#N = 2
k1 = 60
k2 = 40

# Make the betas
beta = matrix(rep(0,num_env*p), nrow = p)
beta[1:k2, ] = 1
if (k1 > 0){
    step_no = 0
    step = k1/num_env
    for (i in 1:num_env){
        beta[k2+1:step+step_no,i] = 1
        step_no = step_no + step
    }
}
beta = beta* amplitude / sqrt(ns[2])
#heatmap(beta,Rowv = NA,Colv = NA, scale = "none")
set.seed(1) ## make sure our model is not changing
permu = sample(1:p)
beta = beta[permu,]
nonnulls = intersect(which(beta[,1] != 0), which(beta[,2] != 0))

fdps = NULL
powers = NULL
fdp = function(selected, nonnulls) length(setdiff(selected, nonnulls)) / max(1, length(selected))
power = function(selected, nonnulls) length(intersect(selected, nonnulls)) / length(nonnulls)                                                   

ptm <- proc.time()
for (iter in 1:N){
    
    Xs = list()
    ys = list()
    for (i in 1:num_env){
        Sigma = toeplitz(rhos[i]^(0:(p-1)))
        X = matrix(rnorm(ns[i]*p),ns[i]) %*% chol(Sigma)
        Xs[[i]] = X
        y = X %*% beta[,i] + rnorm(ns[i])
        ys[[i]] = y
    }
    
    ## create knockoffs
    X_ks = list()
    for (i in 1:num_env){
        X_ks[[i]] = create.gaussian(Xs[[i]], rep(0,p), toeplitz(rhos[i]^(0:(p-1))))
    }
    
    ## W for each env
    W_matrix = NULL
    for (i in 1:num_env){
        W = stat.lasso_coefdiff(Xs[[i]], X_ks[[i]], ys[[i]])
        W_matrix = cbind(W_matrix,W)
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
    thres_stack = knockoff.threshold(W_stack)
    selected_stack = which(W_stack >= thres_stack)
    fdp_stack = fdp(selected_stack, nonnulls)
    power_stack = power(selected_stack,nonnulls)
    
    ## Heuristic
    selected_h = NULL
    for (i in 1:num_env){
        W_i = W_matrix[,i]
        thres_i = knockoff.threshold(W_i, fdr=q)
        selected_h = cbind(selected_h, W_i >= thres_i)
    }
    selected_h = which(as.logical(apply(selected_h, 1, combine_mag, r = 2)))
    fdp_h = fdp(selected_h, nonnulls)
    power_h = power(selected_h,nonnulls)
    
    ## Invariant
    selected_i = invariant_model(W_matrix)
    fdp_i = fdp(selected_i, nonnulls)
    power_i = power(selected_i,nonnulls)
    
    ## Invariant with cross prior statistics
    W_cp_matrix = compute_stats_with_prior(Xs, X_ks, ys, verbose=TRUE)
    selected_c = invariant_model(W_cp_matrix, q=q)
    fdp_c = fdp(selected_c, nonnulls)
    power_c = power(selected_c,nonnulls)
    
    ## Invariant with stack magnitude
    W_sign = sign_fun(apply(W_matrix, 1, min))
    W_eff_w = W_sign * abs(W_stack)
    thres_w = knockoff.threshold(W_eff_w, fdr=q)
    selected_w = which(W_eff_w >= thres_w)
    fdp_w = fdp(selected_w, nonnulls)
    power_w = power(selected_w,nonnulls)
    
    fdps = cbind(fdps, c(fdp_stack, fdp_h, fdp_i, fdp_c, fdp_w))
    powers = cbind(powers, c(power_stack, power_h, power_i, power_c, power_w))
    #print(fdps)
    #print(powers)
    print(iter)
}

fdps
powers


