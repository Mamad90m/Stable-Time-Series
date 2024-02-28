
# Goodness of Fit Tests for Non-Causal Autoregressive Time Series#
# Note: S is the matrix of simulated stable random variable and param1
# and param2 include estimated parameters of AR(1) and AR(2) respectively
# Residuals for non-causal AR(1)#

library(ZIM)
library(portes)
library(WeightedPortTest)

z.hat <- function(x){
  n <- length(x)
  z <- x[1:(n - 1)] - x[n] * bshift(x[1:(n - 1)], k = 1)
  return(z[2:(n - 1)])
}

Z <- t(apply(cbind(S, param1[, 1]), 1, z.hat))

# Residuals for mixed causal and non-causal AR(2)#
 
z.hat <- function(x){
  n <- length(x)
  z <- x[1:(n - 2)] - (x[n - 1] + x[n - 2]) * bshift(x[1:(n - 2)], k = 1) + 
    (x[n - 1] * x[n - 2]) * bshift(x[1:(n - 2)], k = 2) 
  return(z[3:(n - 2)])
}

Z <- t(apply(cbind(S, param2[, c(1, 2)]), 1, z.hat))
#############################################################

## Ljung-Box test (Q_lb)##
Q_lb <- function(x){
  x[x <= sort(x)[length(x) * 0.01]] <- 0
  x[x >= sort(x)[length(x) * 0.99]] <- 0
  
  t <- LjungBox(x)
  p_value <- t[, 4]
  result <-ifelse(p_value > 0.05, 0, 1)
  return(result)
}

T_lb <- t(apply(Z, 1, Q_slb))
rej_lb <- apply(T_lb, 2, sum)   #number of rejected tests
##############################################################

## Monti test (Q_mt)##
Q_mt <- function(x){
  x[x <= sort(x)[length(x) * 0.01]] <- 0
  x[x >= sort(x)[length(x) * 0.99]] <- 0
  
  p_value <- NULL
  for(i in seq(5, 30, 5)) 
    p_value[i] <- as.numeric(Weighted.Box.test(x, lag = i, type = "Monti", 
                                               weighted = FALSE)[[3]])
  p_value <- p_value[!is.na(p_value)]
  result <-ifelse(p_value > 0.05, 0, 1)
  return(result)
}


T_mt <- t(apply(Z, 1, Q_mt))
rej_mt <- apply(T_mt, 2, sum)   #number of rejected tests
###############################################################

## Mahdi-Mc Leod test (Q_mm)##
Q_mm <- function(x){
  x[x <= sort(x)[length(x) * 0.01]] <- 0
  x[x >= sort(x)[length(x) * 0.99]] <- 0
  
  t <- MahdiMcLeod(x)
  p_value <- t[, 4]
  result <-ifelse(p_value > 0.05, 0, 1)
  return(result)
}


T_mm <- t(apply(Z, 1, Q_mm))
rej_mm <- apply(T_mm, 2, sum)   #number of rejected tests
##############################################################

## Fisher-Gallagher test (Q_fg)##
Q_fg <- function(x){
  x[x <= sort(x)[length(x) * 0.01]] <- 0
  x[x >= sort(x)[length(x) * 0.99]] <- 0
  
  p_value <- NULL
  for(i in seq(5, 30, 5)) 
    p_value[i] <- as.numeric(Weighted.Box.test(x, lag = i, type = "Ljung")[[3]])
  p_value <- p_value[!is.na(p_value)]
  result <-ifelse(p_value > 0.05, 0, 1)
  return(result)
}


T_fg <- t(apply(Z, 1, Q_fg))
rej_fg <- apply(T_fg, 2, sum)   #number of rejected tests
#############################################################
## Tests based on squared residuals

## Ljung-Box test (Q_slb)##
Q_slb <- function(x){
  x[x <= sort(x)[length(x) * 0.01]] <- 0
  x[x >= sort(x)[length(x) * 0.99]] <- 0
  
  t <- LjungBox(x, sqrd.res = TRUE)
  p_value <- t[, 4]
  result <-ifelse(p_value > 0.05, 0, 1)
  return(result)
}


T_slb <- t(apply(Z, 1, Q_slb))
rej_slb <- apply(T_slb, 2, sum)   #number of rejected tests
#############################################################

## Monti test (Q_smt)##
Q_smt <- function(x){
  x[x <= sort(x)[length(x) * 0.01]] <- 0
  x[x >= sort(x)[length(x) * 0.99]] <- 0
  
  p_value <- NULL
  for(i in seq(5, 30, 5)){
    p_value[i] <- as.numeric(Weighted.Box.test(x, lag = i, type = "Monti", 
                                               sqrd.res = TRUE, weighted = FALSE)[[3]])
  }
  p_value <- p_value[!is.na(p_value)]
  result <-ifelse(p_value > 0.05, 0, 1)
  return(result)
}


T_smt <- t(apply(Z, 1, Q_smt))
rej_smt <- apply(T_smt, 2, sum)   #number of rejected tests
#############################################################

## Mahdi-Mc Leod test (Q_smm)##
Q_smm <- function(x){
  x[x <= sort(x)[length(x) * 0.01]] <- 0
  x[x >= sort(x)[length(x) * 0.99]] <- 0
  
  t <- MahdiMcLeod(x, sqrd.res = TRUE)
  p_value <- t[, 4]
  result <-ifelse(p_value > 0.05, 0, 1)
  return(result)
}


T_smm <- t(apply(Z, 1, Q_smm))
rej_smm <- apply(T_smm, 2, sum)   #number of rejected tests
#############################################################

## Fisher-Gallagher test (Q_sfg)##
Q_sfg <- function(x){
  x[x <= sort(x)[length(x) * 0.01]] <- 0
  x[x >= sort(x)[length(x) * 0.99]] <- 0
  
  p_value <- NULL
  for(i in seq(5, 30, 5)) 
    p_value[i] <- as.numeric(Weighted.Box.test(x, lag = i, type = "Ljung", 
                                               sqrd.res = TRUE)[[3]])
  p_value <- p_value[!is.na(p_value)]
  result <-ifelse(p_value > 0.05, 0, 1)
  return(result)
}


T_sfg <- t(apply(Z, 1, Q_sfg))
rej_sfg <- apply(T_sfg, 2, sum)   #number of rejected tests
#############################################################










