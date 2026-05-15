#' Sample size required to develop a risk prediction model for binary outcomes
#'
#' @description
#' This function calculates the sample size required to achieve an expected Calibration Slope (S), given anticipated features of the data and the model
#' (outcome prevalence, C-statistic and number of predictors).
#'
#'It takes approximately one minute to run. Ideally it should be followed by checking also
#' the Mean Absolute Prediction Error that corresponds to the calculated sample size.

#' @param S (numeric) The target expected calibration slope
#' @param p (numeric) The anticipated outcome prevalence
#' @param c (numeric) The anticipated C-statistic
#' @param n.predictors (numeric) The number of candidate predictor variables
#' @param nsim (numeric) The number of simulations (at least 500, default value 1000 to ensure small simulation error)
#' @param nval (numeric) Size of validation data (at least 10000 )
#' @param parallel (logical) parallel processing to speed up computations (default=TRUE)
#' @param plot (logical) figures for intermediate plots (default=TRUE)


#'
#' @return n: the required sample size
#'
#' @examples
#' # Find the sample size required for an average calibration slope of S = 0.9
#' # samplesizedev_binary(S = 0.9, p = 0.2, c = 0.85, n.predictors = 10,  nsim = 500, parallel = FALSE)
#'
#' # Prefer parallel computing with >2 cores that ensure faster running
#' # samplesizedev_binary(S = 0.9, p = 0.2, c = 0.85, n.predictors = 10,  nsim = 1000, parallel = TRUE)
#'
#' # Check the expected MAPE and Calibration Slope for the selected size
#' # expected_cs_mape_binary(n = 530, p = 0.2, c = 0.85, n.predictors = 10, nsim = 1000, parallel = TRUE)
#'
#'
#'
#' @seealso
#' expected_cs_mape


samplesizedev_binary_s <- function(S, p, c,   n.predictors, beta = rep(1/n.predictors, n.predictors), nval = 25000, nsim = 1000, parallel = TRUE, plot = TRUE, quick = TRUE, tol){

  set.seed(2022)

  mean_var_eta     <- find_mu_sigma(p, c, tol = 0.00001)
  mean_eta         <- mean_var_eta[1]
  variance_eta     <- mean_var_eta[2]

  # Find beta that corresponds to that variance

  beta    <- beta * sqrt(mean_var_eta[2]/sum(beta^2))
  sigma   <- diag(1, n.predictors)

  # True R2 and R2CS
  MaxR2      <- 1-(((p^(p))*((1-p)^(1-p)))^2)
  ncalc      <- 500000
  x          <- mvtnorm::rmvnorm(ncalc, rep(0, n.predictors), sigma = sigma )
  eta        <- mean_eta+x%*% beta
  y          <- stats::rbinom(ncalc,  1, invlogit(eta))
  p_true     <- as.vector(invlogit(mean_eta + x%*%beta))
  a          <- RcppNumerical::fastLR(cbind(1,x), y)
  L1         <- a$loglikelihood
  L0         <- sum(y*log(mean(y)) + (1-y)*log(1-mean(y)))
  LR         <- -2*(L0-L1)
  r2         <- 1 - exp(-LR/ncalc)


  c_adj  <- adjusted_c_mu_sigma(mean_eta, variance_eta, n.predictors, p, set.seed=1)

  n_init <- round((n.predictors)/ ((S-1)*log(1-  c_adj[2]/S)))

  # r2   <- as.numeric(approximate_R2(c, p, n = 300000)[2])

  n_rvs <- round((n.predictors)/ ((S-1)*log(1-r2/S)))


  if (c<=0.7  )               {inflation_f   <- 1.1  ; min.opt <- n_rvs*0.4}
  if (c>0.7  & c<=0.8 )       {inflation_f   <- 1.5  ; min.opt <- n_rvs*0.7}
  if (c>0.8  & c<=0.85)       {inflation_f   <- 2.1  ; min.opt <- n_rvs*0.8}
  if (c>0.85 & c<=0.9)        {inflation_f   <- 2.8  ; min.opt <- n_rvs*0.9}


  if (c<=0.7  & n.predictors <6)             {inflation_f   <- 0.9 ; min.opt <- n_rvs*0.3}
  if (c>0.7  & c<=0.8  & n.predictors < 8)   {inflation_f   <- 1.8 ; min.opt <- n_rvs*0.4}
  if (c>0.8  & c<=0.85 & n.predictors < 8)   {inflation_f   <- 2.1 ; min.opt <- n_rvs *0.5}
  if (c>0.85 & c<=0.9  & n.predictors < 8)   {inflation_f   <- 2.8 ; min.opt <- n_rvs*0.7}

  max.opt <- inflation_f * n_rvs

  if (plot==FALSE & quick== FALSE) print("Optimisation Started:...") else
    if (plot==TRUE & quick== FALSE)  print("Optimisation Started: check progress on appearing plots...")

  # Automatically adjust number of simulations to ensure MCSE is not too high
  A               <- 2*p*(1-p)*stats::qnorm(c)^2
  app             <- sqrt(1/(A* n_init)+2/(n_init-2) )
  mce             <- round(app/sqrt(nsim),4)

  # if (mce  > tol )   { nsim    <- ceiling(app^2/tol^2/100)*100; print(paste("Note: Number of simulations increased to nsim=",nsim, "to keep Monte Carlo error small", sep="")) }
  # if (mce  > tol & nsim > 3000)   { nsim    <- 3000;
  # #tol <- max(1*mce/1.5,0.0025)}
  # tol <- round(app/sqrt(nsim),4)}

  if (c>=0.85) { min.opt <- 1.1* n_init; max.opt <- 1.3 * n_init} else
  { min.opt <- 0.9 * n_init; max.opt <- 1.05 * n_init}

  s_est <- function(n, nsim=nsim){

    s <-  expected_s_n_binary(n, S = S, mean_eta = mean_eta, variance_eta = variance_eta,  beta = beta, p = p, c = c,
                              n.predictors = n.predictors, nval = nval, nsim = nsim, parallel=parallel, plot = plot, tol = tol)
    s[1] - S
  }


  if (quick ==  FALSE) n   <- bisection(s_est, min.opt, max.opt, tol = tol, nsim = nsim) else n = n_init

  #print(paste("Required sample size: ", n ))

  size               <- NULL
  size$rvs           <- as.vector(n_rvs)
  size$sim           <- as.vector(round(n))

  size

}
