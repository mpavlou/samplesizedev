#########################################
# Find mean and variance of linear predictor for given prevalence
# and C-statistic

invlogit <- function(x) 1/(1+exp(-x))



find_mu_sigma <- function(target.prev, target.c, min.opt = c(-10,0), max.opt = c(0,5), tol = 0.00001){

  pcfun <- function(x){

    #target.prev = 0.2; target.c=0.7; min.opt = c(-7,0.5); max.opt = c(0,14)

    mean      <- x[1]
    variance  <- x[2]

    f1 = function(x) {
      stats::integrate(function(y) {stats::dnorm(x, mean = mean, sd = sqrt(variance)) * stats::dnorm(y, mean = mean, sd = sqrt(variance)) * (1  + exp(-x)) ^ (-1) * (1  + exp(y)) ^ (-1) },
                -Inf, x)$value
    }

    num = stats::integrate(Vectorize(f1), -Inf, Inf)$value

    f2 = function(x) {
      stats::integrate(function(y) {stats::dnorm(x, mean = mean, sd = sqrt(variance)) * stats::dnorm(y, mean = mean, sd = sqrt(variance)) * (1  + exp(-x)) ^ (-1) * (1  + exp(y)) ^ (-1) },
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

  N        <- 2000000
  lp       <- stats::rnorm(N, mean = out[1], sd = sqrt(out[2]))
  p        <- (1 + exp(-lp)) ^ (-1)
  y        <- stats::rbinom(N, 1, prob = p)
  prev     <- mean(y)
  #c        <- as.vector(pROC::roc(y, lp, quiet = TRUE)$auc)
  c        <- quickcstat(y, lp)
  c(out[1], out[2], prev, c)
}

# Check
# round(find_mu_sigma(0.1, 0.65, tol=0.00001),4)

#####################################################





#########################################
# Find scale and variance of linear predictor for given prevalence and C-statistic

#
# find_lambda_sigma <- function(target.prev, target.c, time=5, min.opt = c(0.001,0.02), max.opt = c(0.5,0.5),tol=0.0001){
#
#   pcfunsurv <- function(x){
#
#     #target.prev = 0.2; target.c=0.7; min.opt = c(-7,0.5); max.opt = c(0,14)
#
#     lambda    <- x[1]
#     variance  <- x[2]
#
#     f3     <- function(y) stats::dnorm(y, mean = 0, sd = sqrt(variance)) * (1 - (exp(-time * lambda)) ^ (exp(y)) )
#    #  f3     <- function(y) stats::dnorm(y, mean = 0, sd = sqrt(variance)) * (1  + exp(-y)) ^ (-1)
#
#     prev   <- stats::integrate(f3, - Inf, Inf, subdivisions = 1000L)$value ;prev
#
#     f1 = function(y) {
#       stats::integrate(function(z) {stats::dnorm(y, mean = 0, sd = sqrt(variance)) * stats::dnorm(z, mean = 0, sd = sqrt(variance)) * (1 - (exp(-time * lambda)) ^ exp(y) ) * (exp(-(time-1) * lambda) ^ exp(z)) },
#                        -Inf, y)$value
#     }
#
#     num = stats::integrate(Vectorize(f1), -Inf, Inf)$value; num
#
#     f2 = function(y) {
#       stats::integrate(function(z) {stats::dnorm(y, mean = 0, sd = sqrt(variance)) * stats::dnorm(z, mean = 0, sd = sqrt(variance)) * (1 - (exp(-time * lambda)) ^ exp(y) ) * (exp(-(time-1) * lambda) ^ exp(z)) },
#                        -Inf, Inf)$value
#     }
#
#     denom <- stats::integrate(Vectorize(f2), -Inf, Inf)$value ;denom
#
#     c      <- num/denom ;c
#
#      # c <- NULL
#      # t <- c(1, 2, 5, 10, 20,30, 40, 100, 150, 200, 1000)
#      # for (time in t) {
#      #   c <- c(c, stats::integrate(Vectorize(f1), -Inf, Inf)$value / stats::integrate(Vectorize(f2), -Inf, Inf)$value)
#      # }
#      #
#      # mean(c)
#      #
#      # c <- weighted.mean(c, w = pweibull(t, shape=1, scale=1/lambda) )
#      # c
#
#     abs( c - target.c)^2 + abs(prev - target.prev )^2
#   }
#
#   out      <- stats::optim(par=c(0.01, 0.5), pcfunsurv, c(min.opt, max.opt, tol = tol))$par
#   out
#
#   N        <- 1000000
#   lp       <- stats::rnorm(N, mean = 0, sd = sqrt(out[2]))
#   p        <- (1 - (exp(-time * out[1])) ^ exp(lp) )
#   # y        <- stats::rbinom(N, 1, prob = p)
#   # prev     <- mean(y)
#   # c        <- as.vector(pROC::roc(y, lp, quiet = TRUE)$auc)
#   # c         <-quickcstat(y, lp)
#   c(out[1], out[2], mean(p))
# }
#
# # Check
#
# out <- round(find_lambda_sigma (0.2, 0.7, time=2, tol=0.00001),4)
# out
#
# n      <- 10000000
# thresh <- 0.001
#
# u      <- runif(n)
# lp     <- stats::rnorm(n, mean = 0, sd = sqrt(out[2]))
# t      <- -log(u)/((out[1]) * exp(lp) )
# censor <- rep(1,n)
# censor[t > quantile(t, p = thresh)]=0
# mean(censor)
# quickcindex(censor, t, lp, seed=1)
#
# n      <- 10000
# thresh <- 0.6
#
# u      <- runif(n)
# lp     <- stats::rnorm(n, mean = 0, sd = sqrt(out[2]))
# t      <- -log(u)/((out[1]) * exp(lp) )
# censor <- rep(1,n)
# censor[t > quantile(t, p = thresh)]=0
# mean(censor)
#
# data<-data.frame(t, censor, lp)
#
# a <- cindex(Surv(t, censor) ~ lp, data=data)
# a
#
#
#
#
#
#
#
# thresh <- 0
#
# u      <- runif(n)
# lp     <- stats::rnorm(n, mean = 0, sd = sqrt(out[2]))
# t      <- -log(u)/((out[1]) * exp(lp) )
# censor <- rep(1,n)
# censor[t > quantile(t, p = thresh)]=0
# mean(censor)
#
#
# quickcindex_nocensoring(t,lp)
#
#

 find_sigma <- function(target.c,  min.opt = 0.1, max.opt = 15, tol=0.0001){

 n       <- 50000
 lambda  <- 1

 # Build formula
 measurevar <- "Surv(t, censor)"
 etavar     <- "eta"
 formula    <- stats::as.formula  (paste("Surv(t, censor)", etavar, sep=" ~ "))


   cfun <- function(x) {
     set.seed(2)
     u      <- stats::runif(n)
     eta    <- stats::rnorm(n, mean = 0, sd = sqrt(x))
     t      <- -log(u)/((lambda) * exp(eta) )
     #c      <- 1 - concordance(t ~ eta)$concordance; c
     censor <- rep(1, n)
     data   <- data.frame(t, censor, eta)

     big <- survival::coxph(formula, data = data, x = TRUE, y = TRUE)

     A1  <- pec::cindex(list("Cox X1"= big), formula = formula, data = data, eval.times = 1)

     c <- A1$AppCindex[["Cox X1"]]

     abs(c - target.c)
   }

 out      <- stats::optimize(cfun, lower=min.opt, upper=max.opt, tol = tol)$minimum
 out
 }


find_sigma_quick <- function(target.c,  min.opt = 0.1, max.opt = 8, tol=0.0001){

  n       <- 50000
  lambda  <- 1

  # Build formula
  measurevar <- "survival::Surv(t, censor)"
  etavar     <- "eta"
  formula    <- stats::as.formula(paste( measurevar, etavar, sep=" ~ "))


  cfun <- function(x) {
    set.seed(2)
    u      <- stats::runif(n)
    eta    <- stats::rnorm(n, mean = 0, sd = sqrt(x))
    t      <- -log(u)/((lambda) * exp(eta) )
    c      <- 1 - survival::concordance(t ~ eta)$concordance; c
    abs(c - target.c)
  }

  out      <- stats::optimize(cfun, lower=min.opt, upper=max.opt, tol = tol)$minimum
  out
}


# system.time(x <- find_sigma(target.c=0.9, min.opt=6, max.opt=15))
# x
# n <-50000
# u      <- runif(n)
# eta    <- stats::rnorm(n, mean = 0, sd = sqrt(x))
# t      <- -log(u)/((lambda) * exp(eta) )
# censor <- rep(1,n)
# c      <- 1 - concordance(t ~ eta)$concordance; c
# # Build formula
# measurevar <- "Surv(t, censor)"
# etavar     <- "eta"
# formula    <- as.formula  (paste( measurevar, etavar, sep=" ~ "))
#
# data <- data.frame(t, censor, eta)
#
# #
#  big <- coxph(formula, data=data, x=TRUE, y=TRUE)
#
#  A1  <- pec::cindex(list("Cox X1"= big), formula=Surv(t, censor)~eta, data=data, eval.times=1)
#
#  cest <- A1$AppCindex[["Cox X1"]]
#  cest
#

