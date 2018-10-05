# Generate linear model a la BCH standardized

gen_std_linmod_bch <- function(mu = 1, n = 100, k = 200, alpha_0 = 1, sigma_zeta = 1, sigma_v = 1){
  beta_0 <- c(1/(1:5), rep(0,5), 1/(1:5), rep(0,k-15))
  eta_0 <- c(1/(1:10),rep(0,k-10))
  
  # Generate correlation matrix
  corr_x <- matrix(data=0, nrow = k, ncol = k)
  for(i in 1:k){
    for(j in 1:k){
      corr_x[i,j] <- 0.5**abs(i-j)
    }
  }
  
  x <-rmvnorm(n, mean = rep(0,k), sigma = corr_x)
  zeta <- rnorm(n, mean = 0, sd = sigma_zeta)
  v <- rnorm(n, mean = 0, sd = sigma_v)
  
  d <- x%*%eta_0 + v
  x <- scale(x) # Standardize
  d <- scale(d) # Standardize
  
  y <- mu + alpha_0*d + x%*%beta_0 + zeta
  
  df <- list('d' = d,'X' = data.frame('X' = x), 'Y' = y)
  return(df)
}


