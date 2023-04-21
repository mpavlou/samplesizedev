#########################################
# Find mean and varaince of linear predictor for given prevalnce 
# and C-statistic 

invlogit <- function(x) 1/(1+exp(-x))



find_mu_sigma <- function(target.prev, target.c, min.opt = c(-10,0), max.opt = c(0,5),tol=0.00001){
  
  pcfun <- function(x){
    
    #target.prev = 0.2; target.c=0.7; min.opt = c(-7,0.5); max.opt = c(0,14)
    
    mean      <- x[1]
    variance  <- x[2]
    
    f1 = function(x) {
      stats::integrate(function(y) {stats::dnorm(x, mean = mean, sd = sqrt(variance)) * stats::dnorm(y, mean = mean,sd = sqrt(variance)) * (1  + exp(-x)) ^ (-1) * (1  + exp(y)) ^ (-1) },
                -Inf, x)$value
    }
    
    num = stats::integrate(Vectorize(f1), -Inf, Inf)$value
    
    f2 = function(x) {
      stats::integrate(function(y) {stats::dnorm(x,mean=mean, sd=sqrt(variance)) * stats::dnorm(y, mean = mean,sd = sqrt(variance)) * (1  + exp(-x)) ^ (-1) * (1  + exp(y)) ^ (-1)},
                -Inf, Inf)$value
    }
    
    denom <- stats::integrate(Vectorize(f2), -Inf, Inf)$value
    
    c      <- num/denom
    f3     <- function(x) stats::dnorm(x, mean=mean, sd = sqrt(variance)) * (1  + exp(-x)) ^ (-1)
    
    prev   <- stats::integrate(f3, - Inf, Inf, subdivisions = 1000L)$value
    abs( c - target.c)^2 + abs(prev - target.prev )^2
  }
  
  out      <- stats::optim(par=c(-2.65,0.1), pcfun, c(min.opt, max.opt, tol = tol))$par
  out
  
  N        <- 1000000
  lp       <- stats::rnorm(N, mean = out[1], sd = sqrt(out[2]))
  p        <- (1 + exp(-lp)) ^ (-1)
  y        <- stats::rbinom(N, 1, prob = p)
  prev     <- mean(y)
  #c        <- as.vector(pROC::roc(y, lp, quiet = TRUE)$auc)
  c         <-quickcstat(y, lp)
  c(out[1], out[2], prev, c)
}

# Check
# round(find_mu_sigma(0.174, 0.7, tol=0.00001),4)
