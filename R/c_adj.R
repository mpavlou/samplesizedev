##################################################################

#' Calculate the means and variance of a normal linear predictor given values of C and prevalence

#' @description
#' This function calculates the adjusted C
#'

#' @param target.prev (numeric) The target prevalence
#' @param target.c (numeric) The target aUC
#' @param check (logical) check using a large dataset whether mean and variance correspond to the target values

#' @return   adjusted C

#' @export

#' @examples
#'
#' c_adjusted <- c_adj(target.prev=0.1, target.c=0.85, p=10)
#' c_adjusted
#' samplesizedev(S=0.9, phi = 0.1, c = c_adjusted, p = 10))
#'
#' ##################################################################


c_adj <- function(target.prev, target.c, p=NULL, set.seed=1) {

  n.predictors <- p
  p            <- target.prev
  c            <- target.c

  set.seed(set.seed)
  mean_var <- find_mu_sigma(p, c, tol=0.00001)
  mean     <- mean_var[1]
  variance <- mean_var[2]

  bj       <- sqrt(variance/n.predictors); bj

  n <- 1000000

  x <- mvtnorm::rmvnorm(n, rep(0, n.predictors), sigma = diag(n.predictors) )
  eta  <- cbind(1, x) %*% as.vector(c(mean, rep(bj, n.predictors) ))
  y   <- rbinom(n,1,plogis(eta))

  delta_j <-NULL
  for (i in 1: n.predictors) {
    x0    <- x[,i][y==0]
    x1    <- x[,i][y==1]

    mu0     <- mean(x0)
    mu1     <- mean(x1)

    n0      <- length(x0)
    n1      <- length(x1)

    sigma20 <- var(x0)
    sigma21 <- var(x1)

    sigmac2 <- ( (n0-1)*sigma20 + (n1-1) * sigma21)/ (n0+n1)

    delta_j[i] <- (mu1-mu0)/sigmac2

  }

  # variance_new  <-  sum(delta_j^2)
  if (n.predictors==1)   variance_new  <-  delta_j^2 else
    variance_new  <-  (mean(delta_j))^2*n.predictors; variance_new

  mean_new      <- find_mu(p, variance_new)[1]
  c_new         <- find_c(mean_new, variance_new); c_new

  # ncalc <- 1000000
  # LP    <- stats::rnorm(ncalc, mean_new, sqrt(variance_new))
  # y     <- stats::rbinom(ncalc, 1, plogis(LP))
  #
  # a   <- RcppNumerical::fastLR(cbind(1,LP), y)
  # L1  <- a$loglikelihood
  # L0  <- sum(y*log(mean(y)) + (1-y)*log(1-mean(y)))
  # LR  <- -2*(L0-L1)
  # R2  <- 1 - exp(-LR/ncalc)

  return(round(c_new,3))

}
