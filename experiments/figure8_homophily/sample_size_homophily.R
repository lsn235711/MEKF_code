## Experiments in Figure 8
## Momophily + contagion
## FDR and Power vs sample size
## Figure 8 runs the following with n_sample_size = 100*(1:20) parallelly
n_sample_size = 2000

##################################################
source("MEKF/functions_multienv.R")
source("MEKF/accumulation_test_functions.R")

##------------------------
## Problem parameters
##------------------------
p = 200
num_env = 2
ns = rep(n_sample_size,num_env)
k = 40
k1 = 25
k2 = 25
overlap = 0
wrongs= list((p-k1+1):p, (p - 25 + overlap - k2 + 1):(p - 25 + overlap))
rhos = c(0.5, 0.4)
amplitude = 8

##------------------------
## Make the betas
##------------------------
beta = matrix(rep(0,num_env*p), nrow = p)
nonnulls = 1:k
beta[nonnulls,] = 1
beta = beta* amplitude / sqrt(1024)

##------------------------
## set of nonnulls
##------------------------
## permute the variables
set.seed(1)
permu = sample(1:p)
beta = beta[permu,]
get_new_index = function(v){return (sort(order(permu)[v]))}
nonnulls = get_new_index(nonnulls)
wrongs= lapply(wrongs,FUN = get_new_index)

results = NULL
for (time in 1:100){
    ##------------------------
    ## generate data
    ##------------------------
    Xs = list()
    ys = list()
    group_interests = list()
    for (i in 1:num_env){
        n_temp = ns[i]
        Sigma = toeplitz(rhos[i]^(0:(p-1)))
        X = matrix(rnorm(n_temp*p),n_temp) %*% chol(Sigma)
        y = X %*% beta[,i] + rnorm(n_temp)
        index = which(y > 0)
        group_interest = 2*rbinom(wrongs[[i]], 1, 0.5) - 1
        group_interests[[i]] = group_interest
        X[index, wrongs[[i]]] = t(t(X[index, wrongs[[i]]]) + group_interest)
        Xs[[i]] = X
        ys[[i]] = y
    }
    
    ##------------------------
    ## Compute empirical mean and variance with more samples
    ##------------------------
    Sigma_list = list()
    mu_list = list()
    for (i in 1:num_env){
        n_large = 100000
        Sigma = toeplitz(rhos[i]^(0:(p-1)))
        X = matrix(rnorm(n_large*p),n_large) %*% chol(Sigma)
        y = X %*% beta[,i] + rnorm(n_large)
        index = which(y > 0)
        X[index, wrongs[[i]]] = t(t(X[index, wrongs[[i]]]) + group_interests[[i]])
        Sigma_list[[i]] = var(X)
        mu_list[[i]] = mean(X)
    }
    
    ##------------------------
    ## Create knockoffs
    ##------------------------
    X_ks = list()
    for (i in 1:num_env){
        X_ks[[i]] = create.gaussian(Xs[[i]], mu_list[[i]], Sigma_list[[i]])
    }
    
    ##------------------------
    ## Make W's
    ##------------------------
    W_matrix = NULL
    for (i in 1:num_env){
        W = stat.lasso_coefdiff(Xs[[i]], X_ks[[i]], ys[[i]])
        W_matrix = cbind(W_matrix,W)
    }
    
    fdp = function(selected, nonnulls) length(setdiff(selected, nonnulls)) / max(1, length(selected))
    power = function(selected, nonnulls) length(intersect(selected, nonnulls)) / length(nonnulls)                                                   
    ##------------------------
    ## MEKF
    ##------------------------
    selected_seq = partial_conjunction(W_matrix, r = 2, method = "seqstep", c = 0.5)
    selected_accu = partial_conjunction(W_matrix, r = 2, method = "accumulation", randomness = 3)
    fdp_accu = fdp(selected_accu, nonnulls)
    power_accu = power(selected_accu, nonnulls)
    fdp_seq = fdp(selected_seq, nonnulls)
    power_seq = power(selected_seq, nonnulls)
    
    ##------------------------
    ## Intersection
    ##------------------------
    selectedi = 1:p
    for (i in 1:num_env){
        thres = knockoff.threshold(W_matrix[,i])
        selectedi1 = which(W_matrix[,i] >= thres)
        selectedi = intersect(selectedi, selectedi1)
    }
    fdpi = fdp(selectedi, nonnulls)
    poweri = power(selectedi, nonnulls)
    
    ##------------------------
    ## Pool the environments
    ##------------------------
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
    
    
    
    ##------------------------
    ## InvariantCausalPrediction
    ##------------------------
    # library(InvariantCausalPrediction)
    # icp = ICP(X_stack, y_stack, ExpInd = rep(1:num_env, time = ns), selection = "lasso")
    # selected_icp = Reduce(union,icp$acceptedSets)
    # fdpicp = fdp(selected_icp, nonnulls)
    # powericp = power(selected_icp,nonnulls)
    result_new = data.frame(fdps = c(fdp_seq, fdp_accu, fdpi, fdp_stack),
                            powers = c(power_seq, power_accu, poweri, power_stack),
                            method = c("Seqstep", "accumulation", "intersection", "stack"),
                            n_sample_size = rep(n_sample_size,4))
    results = rbind(results, result_new)
}

results
