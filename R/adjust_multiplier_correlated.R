
# Find coefficients  for desired c-statistic and prevalence (binary predictors)

adjust_multiplier_correlated <- function(c, mean, beta, n.predictors, min.opt = variance/3, max.opt = variance*3, tol=0.01, cor0, cor1){

  N <- 3000000

  n.noise <- length(beta[beta==0])
  n.true  <- n.predictors - n.noise

  # mean_var <- find_mu_sigma(p, c, tol=0.001)

  # beta_new <- find_multiplier(sqrt(variance), beta=beta, n.true=n.true, n.noise=n.noise)

  # Specify correlation matrix
  sigma <- matrix(0, nrow = n.predictors,  ncol = n.predictors)
  sigma[1:n.true, 1:n.true] <- cor0
  if (n.noise>0) {
    sigma[(n.true+1):n.predictors, (n.true+1):n.predictors] <- cor1}
  diag(sigma) <- 1

  x     <- mvtnorm::rmvnorm(round(N), rep(0,nrow(sigma)), sigma = sigma )

  cfun <- function(adjust){

    eta   <- mean + x %*% as.matrix(beta) * adjust
    y     <- stats::rbinom(N,  1, invlogit(eta))
    cest  <- quickcstat(y,eta)

    abs(cest- c)
  }
  adjust <- optimize(cfun, c(0, 7, tol = 0.0005))$minimum
  adjust

  beta_new <- beta*adjust
  beta_new

   # Check beta's are ok
   #eta <- mean + x%*%as.matrix(beta_new)
   #y     <- rbinom(N,  1, invlogit(eta))
   #mean(y)
   #quickcstat(y, eta)

}

#####################################################################

# Find coefficients  for desired C and prevalnce (binary predictors)

