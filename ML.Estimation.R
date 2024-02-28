
#Maximum Likelihood Estimation for Non-Causal Alph-Stable Autoregressive Models#
### Non-Causal AR(1) ###


library(stabledist)
library(ZIM)
library(MARX)


alpha <- 1.5
beta <- 0     
gamma <- 1
mu <- 0

dist.eps <- c('stable', alpha, betta, gamma, mu) 
obs <- 500
c_par <- 0
nc_par <- 0.5 # AR parameter phi = 2 in this case
n <- 100  # number of realizations

S <- matrix(rep(NA, n * obs), n, obs)
for(i in 1:nrow(S))
  S[i, ] <- sim.marx(dist.eps, NULL, obs, c_par, nc_par, NULL)$y

param1 <- matrix(rep(NA, 5 * n), n, 5)


for(i in 1:nrow(param1))
  param1[i, ] <- as.numeric(unlist(optim(c(1.2, 0.5, 0, 1, 0), function(par){
    ifelse(par[2]<0.1 | par[2]>2 | par[3]< -1 | par[3]>1 | par[4]<=0, 
           -10000, sum(log(dstable((S[i, ] - par[1] * bshift(S[i, ], k = 1))[2:obs], par[2], 
                                   par[3], par[4], par[5]))+ log(abs(par[1]))))
  }, control = list(fnscale = -1, maxit = 2000))))[1:5]


### Mixed Causal and Non-Causal AR(2)###

alpha <- 1.5
beta <- 0     
gamma <- 1
mu <- 0

dist.eps <- c('stable', alpha, beta, gamma, mu) 
obs <- 500
c_par <- 0.5   #phi_1 = 0.5
nc_par <- 0.5  #phi_2 = 2


S <- matrix(rep(NA, n * obs), n, obs)
for(i in 1:nrow(S))
  S[i, ] <- sim.marx(dist.eps, NULL, obs, c_par, nc_par, NULL)$y



param2 <- matrix(rep(NA, 6 * n), n, 6)


for(i in 1:nrow(param2))
  param2[i, ] <- as.numeric(unlist(optim(c(0.5, 1.2, 0.5, 0, 10, 10), function(par){
    ifelse(par[3]<0.1 | par[3]>2 | par[4]< -1 | par[4]>1 | par[5]<=0, 
           -10000, sum(log(dstable(((S[i, ] - (par[1] + par[2])* bshift(S[i, ], 
                                                                        k = 1) + par[1] * par[2] * bshift(S[i, ], k = 2))[3:obs]), par[3], 
                                                par[4], par[5], par[6]))+ log(abs(par[2]))))
  }, control = list(fnscale = -1, maxit = 2000))))[1:6]


