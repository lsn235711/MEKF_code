## Experiments in Figure A1 (first row)
## Consistency and CRT
## Consistency testing
## Should run the following parallelly (e.g. pick a smaller N) to save computation time

##################################################
library(glmnet)

# Problem parameters
q = 0.1
num_env = 3
ns = rep(200, num_env)
p = 100
rhos = rep(0.5,num_env)
amplitude = 3
k1 = 20
k2 = 20
k3 = 20
N = 1000
B = 100

## create beta
beta = matrix(rep(0,num_env*p), nrow = p)
beta[1:k3,] = 1
if (k2 > 0){
    beta[k3 + 1:k2,1:2] = 1
}
if (k1 > 0){
    beta[k3 + k2 + 1:k1,1] = 1
}
beta = beta* amplitude / sqrt(ns[1])


## perturb data
set.seed(1) ## make sure our model is not changing
permu = sample(1:p)
beta = beta[permu,]
nonnulls = which(beta[,1]*beta[,2]*beta[,3] != 0)
## pick a index of interest
indices = order(permu)[c(1, k3+1, k2+k3+1, k1+k2+k3+1)]
no_nonnull_value = c(3,2,1,0)

start_time <- Sys.time()
pjs = NULL
no_nonnulls = NULL
seeds = NULL
for (time in 1:N){
    ## generate data
    Xs = list()
    ys = list()
    for (i in 1:num_env){
        Sigma = toeplitz(rhos[i]^(0:(p-1)))
        X = matrix(rnorm(ns[i]*p),ns[i]) %*% chol(Sigma)
        Xs[[i]] = X
        y = X %*% beta[,i] + rnorm(ns[i])
        ys[[i]] = y
    }
    
    for (jj in 1:4){
        j = indices[jj]
        no_nonnull = no_nonnull_value[jj]
        ps = NULL
        for (env in 1:num_env){
            Sigma = toeplitz(rhos[i]^(0:(p-1)))
            mu = rep(0,p)
            
            ## compute T0
            X = Xs[[env]]
            y = ys[[env]]
            fit = cv.glmnet(X, y)
            lambda = fit$lambda.min
            T0 = abs(coef(fit,s=lambda)[j+1])
            
            ## sample from conditional distribution
            Sigma = toeplitz(rhos[i]^(0:(p-1)))
            Sigma22 = Sigma[-j,-j]
            Sigma22Inv = solve(Sigma22)
            Sigma12 = Sigma[j,-j]
            Sigmacho = sqrt(as.vector((1 - t(Sigma12) %*% Sigma22Inv %*% Sigma12)))
            SigmaProd = Sigma12 %*% Sigma22Inv
            mu_c = (X[,-j] - mu[-j]) %*%  t(SigmaProd) + mu[j]
            
            Ts = NULL
            for (b in 1:B){
                Xnew = X
                Xnew[,j] = Sigmacho * rnorm(ns[i]) + mu_c
                fit = cv.glmnet(Xnew, y)
                lambda = fit$lambda.min
                Ts = c(Ts, abs(coef(fit,s=lambda)[j+1]))
            }
            n_eq = sum(c(T0, Ts) == T0)
            n_large = sum(c(T0, Ts) > T0)
            if(n_eq > 0){
                n_eq = sample(1:n_eq,1)
            }
            ps = c(ps, (n_eq + n_large)/B)
        }
        p1 = max(ps)
        pjs = c(pjs, p1)
        no_nonnulls = c(no_nonnulls, no_nonnull)
    }
}

results

