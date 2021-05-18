
## Experiments in Figure 7
## Biased sampling
## FDR and Power vs overlap
## Figure 7 runs the following with overlap = 0,5,10,15,20,25 parallelly

overlap = 10
##################################################
source("MEKF/functions_multienv.R")
source("MEKF/accumulation_test_functions.R")

##------------------------
## Problem parameters
##------------------------
p = 200
num_env = 2
ns = rep(1200,num_env)
k = 25
k1 = 25
k2 = 25
wrongs= list((p-k1+1):p, (p - 25 + overlap - k2 + 1):(p - 25 + overlap))
rhos = c(0.5, 0.4)
amplitude = 4

##------------------------
## Make the betas
##------------------------
beta = matrix(rep(0,num_env*p), nrow = p)
nonnulls = 1:k
beta[nonnulls,] = 1
beta = beta* amplitude / sqrt(900)

##------------------------
## set of nonnulls
##------------------------
set.seed(1)
## permute the variables
permu = sample(1:p)
beta = beta[permu,]
get_new_index = function(v){return (sort(order(permu)[v]))}
nonnulls = get_new_index(nonnulls)
wrongs= lapply(wrongs,FUN = get_new_index)

fdp = function(selected, nonnulls) length(setdiff(selected, nonnulls)) / max(1, length(selected))
power = function(selected, nonnulls) length(intersect(selected, nonnulls)) / length(nonnulls)   

results = NULL
for (time in 1:1000){
    ##------------------------
    ## generate data
    ##------------------------
    Xs = list()
    ys = list()
    for (i in 1:num_env){
        n_temp = 3*ns[i]
        Sigma = toeplitz(rhos[i]^(0:(p-1)))
        X = matrix(rnorm(n_temp*p),n_temp) %*% chol(Sigma)
        y = X %*% beta[,i] + rnorm(n_temp)
        index = which(y/sqrt(amplitude^2*k/400+1)/2 + rowSums(X[,wrongs[[i]]])/sqrt(k1) >= 0)
        index_sub = index[1:ns[i]]
        X_sub = X[index_sub,]
        y_sub =  y[index_sub]
        Xs[[i]] = X_sub
        ys[[i]] = y_sub
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
        index = which(y/sqrt(amplitude^2*k/400+1)/2 + rowSums(X[,wrongs[[i]]])/sqrt(k1) >= 0)
        X_sub = X[index,]
        Sigma_list[[i]] = var(X_sub)
        mu_list[[i]] = mean(X_sub)
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
                            overlap = rep(overlap,4))
    results = rbind(results, result_new)
}

results
