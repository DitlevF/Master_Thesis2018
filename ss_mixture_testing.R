# Test mixture algorithm

rm(list=ls())
library(coda)

# Set Working Directory

setwd('/Users/ditlevkf/OneDrive/R Codes/')

# Load programs
source('gen_std_linmod_bch.R')
source('ss_mixture_vol1.R')

# Load data II: a la BCH
nobs <- 100
kvar <- 50 # Must be larger than 25

df <- gen_std_linmod_bch(n = nobs, k = kvar, alpha_0 = 1, sigma_zeta = 1, sigma_v = 1 )
d <- df$d; x <- df$X; Y <- df$Y
X <- as.matrix(cbind(d,x))

beta_true <- c(rep(1,6),rep(0,5),rep(1,5),rep(0,kvar-15)) # TRUE/FALSE vector of regressors.

draws <- ss_mixture(X,Y, iter = 1000, c = 10, fix = c(0,0,0,1, rep(0,kvar-3))); # Note, higher c, higher exclusion probability.

alphas <- draws$alphas

deltas <- draws$deltas
mu <- draws$mu
sigma2 <- draws$sigma2
omega <- draws$omega
z <- draws$z

# Evaluate MCMC
post_alpha <- as.mcmc(alphas)
plot(post_alpha)

delta_mat <- as.matrix(deltas)
colMeans(delta_mat)
