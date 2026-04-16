adjusted_c_mu_sigma <-function(mean, variance, n.predictors, p, set.seed=1) {

  set.seed(set.seed)
  # mean_var <- find_mu_sigma(p, c, tol=0.00001)
  # mean     <- mean_var[1]
  # variance <- mean_var[2]

  bj       <- sqrt(variance/n.predictors); bj
  # cj       <- find_c(mean, bj^2) ;cj
  # mean_varj <- find_mu_sigma(p, cj, tol=0.00001)
  #
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

  c_new        <- find_c(mean_new, variance_new); c_new

  ncalc <- 1000000
  LP    <- stats::rnorm(ncalc, mean_new, sqrt(variance_new))
  y     <- stats::rbinom(ncalc, 1, plogis(LP))

  # Fit a logistic regression with LP as covariate;
  # this is essentially a calibration model, and the intercept and
  # slope estimate will ensure the outcome proportion is accounted
  # for, without changing C-statistic

  # fit <- rms::lrm(y~LP)
  #
  # max_R2 <- function(prev){
  #   1-(prev^prev*(1-prev)^(1-prev))^2
  # }
  #
  # R2 <- as.numeric(fit$stats['R2']) * max_R2(p)

  a   <- RcppNumerical::fastLR(cbind(1,LP), y)
  L1  <- a$loglikelihood
  L0  <- sum(y*log(mean(y)) + (1-y)*log(1-mean(y)))
  LR  <- -2*(L0-L1)
  R2  <- 1 - exp(-LR/ncalc)

  c(c_new, R2)

}

find_mu <- function(target.prev, variance, min.opt = c(-10), max.opt = c(0.02), tol = 0.0001){


  pfun <- function(x){

    mean      <- x[1]
    f3     <- function(x) stats::dnorm(x, mean=mean, sd = sqrt(variance)) * (1  + exp(-x)) ^ (-1)
    prev   <- stats::integrate(f3, - Inf, Inf, subdivisions = 1000L)$value

    abs(prev - target.prev)
  }

  # out      <- stats::optim(par=c(-1), pfun, c(min.opt, max.opt, tol = tol))$par

  out      <- stats::optimize(pfun, c(min.opt, max.opt), tol = tol)$min
  out


  N        <- 500000

  #Better 2000000 to check
  lp       <- stats::rnorm(N, mean = out[1], sd = sqrt(variance))
  p        <- (1 + exp(-lp)) ^ (-1)
  y        <- stats::rbinom(N, 1, prob = p)
  prev     <- mean(y)
  c        <- quickcstat(y, lp)
  c(out[1], variance, prev, c)
}



adjusted_c<-function(mean, variance, n.predictors=NULL, set.seed=1) {

  set.seed(set.seed)
  mean_var <- find_mu_sigma(p, c, tol=0.00001)
  mean     <- mean_var[1]
  variance <- mean_var[2]

  bj       <- sqrt(variance/n.predictors); bj
  # cj       <- find_c(mean, bj^2) ;cj
  # mean_varj <- find_mu_sigma(p, cj, tol=0.00001)
  #
  n <- 2000000

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

  # c_new        <- find_c(mean_new, variance_new); c_new

  ncalc <- 1000000
  LP    <- stats::rnorm(ncalc, mean_new, sqrt(variance_new))
  y     <- stats::rbinom(ncalc, 1, plogis(LP))

  # Fit a logistic regression with LP as covariate;
  # this is essentially a calibration model, and the intercept and
  # slope estimate will ensure the outcome proportion is accounted
  # for, without changing C-statistic

  # fit <- rms::lrm(y~LP)
  #
  # max_R2 <- function(prev){
  #   1-(prev^prev*(1-prev)^(1-prev))^2
  # }
  #
  # R2 <- as.numeric(fit$stats['R2']) * max_R2(p)
  #

  a   <- RcppNumerical::fastLR(cbind(1,LP), y)
  L1  <- a$loglikelihood
  L0  <- sum(y*log(mean(y)) + (1-y)*log(1-mean(y)))
  LR  <- -2*(L0-L1)
  R2  <- 1 - exp(-LR/ncalc)

  c(c_new,R2)

}


find_c <- function(mean, variance) {

  f1 = function(x) {
    stats::integrate(function(y) {stats::dnorm(x, mean = mean, sd = sqrt(variance)) * stats::dnorm(y, mean = mean, sd = sqrt(variance)) * (1  + exp(-x)) ^ (-1) * (1  + exp(y)) ^ (-1) },
                     -Inf, x)$value
  }


  f2 = function(x) {
    stats::integrate(function(y) {stats::dnorm(x, mean = mean, sd = sqrt(variance)) * stats::dnorm(y, mean = mean, sd = sqrt(variance)) * (1  + exp(-x)) ^ (-1) * (1  + exp(y)) ^ (-1) },
                     -Inf, Inf)$value
  }

  f3     <- function(x) stats::dnorm(x, mean=mean, sd = sqrt(variance)) * (1  + exp(-x)) ^ (-1)


  num = stats::integrate(Vectorize(f1), -Inf, Inf)$value

  denom <- stats::integrate(Vectorize(f2), -Inf, Inf)$value

  stats::integrate(Vectorize(f3), -Inf, Inf)$value

  cest     <- num/denom;

  cest

}


invlogit <- function(x) 1/(1+exp(-x))

logit <-function(x) log(x/(1-x))



