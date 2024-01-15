
#' Sample size required to develop a risk prediction model for Survival outcomes
#'
#' @description
#' This function calculates the sample size required to achieve an expected Calibration Slope (S), given anticipated features of the data and the model
#' (proportion of events (uncensored observations) at a given time point of interest, C-index and number of predictors).
#'
#'It takes approximately one minute to run. Ideally it should be followed by checking also
#' the Mean Absolute Prediction Error that corresponds to the calculated sample size.

#' @param S (numeric) The target expected calibration slope
#' @param p (numeric) The anticipated proportion of events (at the time-point of interest)
#' @param c (numeric) The anticipated C-index
#' @param n.predictors (numeric) The number of candidate predictor variables
#' @param beta (numeric) The relative strength of predictors (0 for noise)
#' @param nsim (numeric) The number of simulations (at least 500, default value 1000 to ensure small simulation error)
#' @param nval (numeric) Size of validation data (at least 10000 )
#' @param parallel (logical) parallel processing to speed up computations (default=TRUE)

#'
#' @return n: the required sample size
#'
#' @examples
#' # Find the sample size
#' # samplesizedev_survival(S = 0.9, p = 0.2, c = 0.85, n.predictors = 10,  nsim = 500, parallel = FALSE)
#'
#' # Prefer parallel computing with >2 cores that ensure faster running
#' # samplesizedev_survival(S = 0.9, p = 0.2, c = 0.85, n.predictors = 10,  nsim = 1000, parallel = TRUE)
#'
#' # Check the expected MAPE and Calibration Slope for the selected size
#' # expected_cs_survival(n = 530, p = 0.2, c = 0.85, n.predictors = 10, nsim = 1000, parallel = TRUE)
#'
#'
#'
#' @seealso
#' expected_cs_mape_binary


samplesizedev_survival <- function(S, p, c,  n.predictors, beta, nval = 25000, nsim = 1000, parallel = TRUE){

  set.seed(1)
  lambda <-1

  if (c<=0.8) variance_eta <-  find_sigma_quick(c)[1]  else variance_eta <-  find_sigma(c)[1]
  p.censor <- 1 - p

  # check
  # Size - big data

  #beta      <- rep(1, n.predictors)
  f         <- sqrt(variance_eta / sum(beta^2))
  beta      <- f * beta

  sigma     <- diag(1, length(beta))
  x         <- mvtnorm::rmvnorm(nval, rep(0, length(beta)), sigma = sigma)
  eta       <- x %*% as.vector(beta)


  u         <- stats::runif(nval)
  t         <- -log(u)/( (lambda) * exp(eta) )

  censor   <- rep(1,nval)
  cutoff   <- stats::quantile(t, 1- p.censor)
  censor[t > cutoff]=0
  ptrue    <- mean(censor) ;ptrue; cutoff

  data <- data.frame(t, censor, eta)

  # Build formula
  xvars      <- paste( "X", seq(1 : n.predictors), sep = "")
  measurevar <- "survival::Surv(t, censor)"
  etavar     <- "eta"
  formula    <- stats::as.formula  (paste( measurevar, paste(xvars, collapse=" + "), sep=" ~ "))
  formula    <- stats::as.formula  (paste( measurevar, etavar, sep=" ~ "))

  big <- survival::coxph(formula, data=data, x=TRUE, y=TRUE)

  #A1  <- pec::cindex(list("Cox X1"= big),
  ##                  formula=formula,
  #                  data=data,
  #                  eval.times=2)

  #cest <- A1$AppCindex[["Cox X1"]]
  #cest


  logtest   <- -2 * (big$loglik[1] - big$loglik[2])
  rsq       <- c(rsq = 1 - exp(-logtest/big$n), maxrsq = 1 -
                   exp(2 * big$loglik[1]/big$n))
  r2cens    <- rsq[1]
  r2cens
  mean(censor)


  n_init  <- ceiling((n.predictors)/ ((S-1)*log(1-r2cens/S))); n_init


  min.opt                              <- round(n_init*0.8)
  if (c<=0.7)            inflation_f   <- 1.5
  if (c>0.7  & c<=0.8)   inflation_f   <- 2
  if (c>0.8  & c<=0.85)  inflation_f   <- 3.5
  if (c>0.85 & c<=0.9)   inflation_f   <- 4
  max.opt                              <- round(inflation_f*n_init)

  tol = ceiling(round(n_init/100)/10)*10
  #tol = 20

  print("Optimisation Starting ~ 4 min left...")

  s_est <- function(n, nsim = nsim) {

    s <-  expected_s_n_survival(n, S = S, variance_eta = variance_eta,  p = p, c = c, beta = beta, n.predictors = n.predictors, nval = nval, nsim = nsim, parallel=parallel)
    s[1] - S
  }

  n <- bisection(s_est, min.opt, max.opt, tol = tol, nsim = nsim)
  tol = ceiling(round(n_init/200)/5) * 5
  n <- ceiling(n/tol)*tol

  #run <- expected_s(n, p=p, c=c, n.true=n.true, n.noise=n.noise, beta = c(0.5,0.3,0.3,0.15,0.15), nsim=1000, nval=50000, cores=2)

  size        <- NULL
  size$rvs  <- as.vector(n_init)
  size$sim <- as.vector(n)

  size

}


# system.time ( as <- samplesizedev_survival(S = 0.9, p = 0.5, c = 0.75, n.predictors = 12,  nsim = 1000, parallel = TRUE, nval = 10000))
#
#
 #system.time(e_actual <- expected_cs_survival(n = as$actual, p = 0.5, c = 0.75, n.predictors = 12, nsim = 1000, nval = 25000, parallel = TRUE))
# e_actual
#
# system.time(e_riley <- expected_cs_survival(n = as$riley, p = 0.5, c = 0.75, n.predictors = 12, nsim = 1000, nval = 25000, parallel = TRUE))
# e_riley
#
#
# system.time ( ab<- samplesizedev(S = 0.9, p = 0.1, c = 0.85, n.predictors = 12,  nsim = 1000, parallel = TRUE, nval = 25000))
# n <-as.numeric(gsub("\\D", "", ab))
#
# system.time ( eab<- expected_cs_mape(n = n, p = 0.1, c = 0.85, n.predictors = 12, nsim = 1000, nval = 25000, parallel = TRUE))
# eab






