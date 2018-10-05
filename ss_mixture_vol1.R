
# Back-up
ss_mixture <- function(X, Y, a =  c(1,2,1,8,5,7,8,9,10,12,14,1000000000), 
                        b = c(3,3,2,8,1,1,1,1,1,1,1,1), fix = rep(0,ncol(X)), c = 1, iter = 1000, print_iter = TRUE){
  require(mvtnorm)
  require(MCMCpack)
  
  groups <- length(a)
  k <- length(X[1,])
  N <- length(X[,1])
  
  # Priors for independence slab
  a0 <- rep(0,k)
  A0 <- diag(k)*c
  
  # Matrices to store chains
  delta_draws <- matrix( data = NA, nrow = iter, ncol = k)
  sigma_draws <- rep(NA, iter)
  mu_draws <- rep(NA, iter)
  alpha_draws <- matrix(data = NA, nrow = iter, ncol = k)
  omega_draws <- matrix(data = NA, nrow = iter, ncol = k)
  z_draws <- matrix(data = NA, nrow = iter, ncol = k) # Groups draws
  
  # Initialize delta, omega and alpha (and gamma and p for sampling omega)
  delta <- rep(TRUE,k) == fix # Initialize to include deltas we want in the model, rest out.
  alpha <- rep(0,k)
  
  gamma <- rep(1,groups) # There are two groups
  omega <- runif(k); omega[fix == 1] <- 1 # Initialize omega
  p <- rdirichlet(1,gamma) # Prior inclusion probability initialization
  
  # Test algorithm
  Rj_draws <- matrix(data = NA, nrow = iter, ncol = k)
  post_delta_draws <- matrix(data = NA, nrow = iter, ncol = k)
  
  for(j in 1:iter){ 
    
    if(print_iter == TRUE){# Print every 100 iteration 
      if (j %%100 == 0){ 
        cat("iter =", j, "\n")
      }
    }
    # Step Ia: Deltas
    perm <- sample(1:k,k) # Draw a random permutation
    for(i in 1:k){
      log_yd <- rep(NA,2)
      log_delta <- rep(NA,2)
      log_post_delta <- rep(NA,2)
      
      for(d in 0:1){ # Compute likelihoods for delta_j = 0 and delta_j = 1 given delta_\j
        delta[perm[i]] <- d == 1 # Alternate between TRUE/FALSE
        log_yd[d+1] <- ml_y(a0,A0,Y,X,delta)
      }
      
      R_j <- exp(log_yd[1] - max(log_yd))/exp(log_yd[2] - max(log_yd)) # Subtract max(log_yd) to avoid very low numbers
      post_delta <- c(1-1/(1+R_j*((1-omega[perm[i]])/omega[perm[i]])),1/(1+R_j*((1-omega[perm[i]])/omega[perm[i]])))
      delta[perm[i]] <- sample(x=c(0,1),1, replace = TRUE, prob = post_delta)
      delta <- delta == TRUE # Convert to TRUE/FALSE
      
      # Test algorithm
      Rj_draws[j,perm[i]] <- R_j
      post_delta_draws[j, perm[i]] <- post_delta[1]
    }
    delta_draws[j,] <- delta
    
    # Auxiliary step: get aN, AN, sN, SN
    c <- get_aN_AN_sn_SN(a0,A0,Y,X,delta)
    aN <- c$aN; AN <- c$AN; sN <- c$sN; SN <- c$SN
    
    # Step Ib: sigmas
    
    sigma2 <- 1/rgamma(1,sN, SN)
    sigma_draws[j] <- sigma2
    
    # Step Ic: mu
    mu <- rnorm(1, mean = mean(Y), sd = sqrt(sigma2/N))
    mu_draws[j] <- mu
    
    # Step II: omega
    df <- sample_omega(delta, omega, p, a = a,  b = b)
    omega <- df$omega; p <- df$p; z <- df$z
    omega_draws[j,] <- omega 
    z_draws[j,] <- z
    
    # STEP III: Non-zero alphas
    alpha[delta] <- rmvnorm(1, mean = aN, sigma = sigma2 * AN )
    alpha[delta == FALSE] <- 0
    alpha_draws[j,] <- alpha
  }
  
  draws <- list('alphas' = data.frame('alphas' = alpha_draws), 'deltas' = data.frame('deltas' = delta_draws), 
                'mu' = mu_draws, 'sigma2' = sigma_draws, 
                'omegas' = data.frame('omegas' = omega_draws), 'Rj' = data.frame('Rj' = Rj_draws), 
                'post_delta' = data.frame('post_delta' = post_delta_draws), 'z' = data.frame('z' = z_draws))
  
  colnames(draws$deltas) <- colnames(draws$alphas) <- colnames(X)
  
  return(draws)
  
}

#######################################################################
#######################################################################
####################### Auxiliary functions ###########################
#######################################################################
#######################################################################

ml_y <- function(a0, A0, y, X, delta){
  N <- length(X[,1])
  yc <- y-mean(y)
  sN <- (N-1)/2
  
  if(sum(delta)>0){
    X_d <- X[,delta]
    A0_d <- as.matrix(A0[delta,delta]) # Make sure it's matrix format, if it's e.g. a scalar
    a0_d <- a0[delta]
    
    inv_A0 <- solve(A0_d)
    AN_d <- solve(t(X_d) %*% X_d + inv_A0)
    aN_d <- AN_d %*% (t(X_d) %*% yc + inv_A0 %*% a0_d)
    SN_d <- 0.5 * (t(yc) %*% yc + t(a0_d) %*% inv_A0 %*% a0_d - t(aN_d) %*% solve(AN_d) %*% aN_d)
    
    log_ml_yd <- -0.5*log(N)-0.5*(N-1)*(log(2*pi))+0.5*(log(det(AN_d)))-0.5*(log(det(A0_d)))+ lgamma(sN)-sN*(log(SN_d))
  }
  else{
    SN_d <- 0.5*(t(yc)%*%yc)
    log_ml_yd <- -0.5*log(N)-0.5*(N-1)*(log(2*pi)) + lgamma(sN) - sN * (log(SN_d))
  }
  return(log_ml_yd)
}

get_aN_AN_sn_SN <- function(a0,A0,y,X,delta){
  N=length(X[,1])
  yc <- y-mean(y)
  sN <- (N-1)/2
  
  if(sum(delta)>0){
    X_d <- X[,delta]
    A0_d <- as.matrix(A0[delta,delta]) # Make sure it's matrix format
    a0_d <- a0[delta]
    inv_A0 <- solve(A0_d)
    
    AN_d <- solve(t(X_d) %*% X_d + inv_A0)
    aN_d <- AN_d %*% (t(X_d) %*% yc + inv_A0 %*% a0_d)
    SN_d <- 0.5 * (t(yc) %*% yc + t(a0_d) %*% inv_A0 %*% a0_d - t(aN_d) %*% solve(AN_d) %*% aN_d)
  }
  else{
    AN_d <- matrix(1)
    aN_d <- c(1)
    SN_d <-0.5*(t(yc)%*%yc)
  }
  return(list('aN' = aN_d, 'AN' = AN_d, 'sN' = sN, 'SN' = SN_d))
}


# Sample Omega

sample_omega <- function(delta, omega, p, a = c(1,1,1,8,5,7,8,9,10,12,14,1000000000), b = c(4,3,2,8,1,1,1,1,1,1,1,1)){
  require(MCMCpack)
  
  k <- length(delta)
  groups <- length(a)
  
  # Prior on groups
  gamma <- rep(1, groups)
  
  # Sample z
  
  prob <- matrix(NA, nrow = groups, ncol = k) # Each column j represents omega_j's probability of belonging to group i (row i)
  z <- matrix(NA, nrow = groups, ncol = k) 
  
  for(i in 1:groups){
    prob[i,] <- p[i] * dbeta(omega, a[i], b[i]) # omega^(a[i]-1)*(1-omega)^(b[i]-1) / beta(a[i],b[i])
  }
  #print(prob)
  for(j in 1:k){
    prob[,j] <- prob[,j] / sum(prob[,j])
    #print(prob)
    z[,j] <- as.integer(rmultinom(1,1,prob[,j]))
  }
  #print(prob)
  #print(z)
  n <- rep(NA,groups)
  for(i in 1:groups) n[i] <- sum(z[i,])
  
  Z <- rep(NA,k) # Create categorical group category, e.g. Z[j] = 2 when omega_j belongs to group 2
  for(j in 1:k) Z[j] <- which(z[,j] == 1)
  
  # Sample p
  p <- rdirichlet(1, gamma + n)
  
  # Sample omega
  
  for(j in 1:k){
    omega[j] <- rbeta(1, a[Z[j]] + delta[j], b[Z[j]] + 1 - delta[j])
  }
  
  return(list("omega" = omega, "p" = p, "z" = Z))
}


