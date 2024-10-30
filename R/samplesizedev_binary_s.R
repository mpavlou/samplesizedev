
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


samplesizedev_binary_s <- function(S, p, c,   n.predictors, beta = rep(1/n.predictors, n.predictors), nval = 25000, nsim = 1000, parallel = TRUE){

  set.seed(2022)

  mean_var_eta     <- find_mu_sigma(p, c, tol = 0.00001)
  mean_eta         <- mean_var_eta[1]
  variance_eta     <- mean_var_eta[2]

  r2   <- as.numeric(approximate_R2(c, p, n = 300000)[2])

  n_init <- round((n.predictors)/ ((S-1)*log(1-r2/S)))


  if (c<=0.7  )               {inflation_f   <- 1.1 ; min.opt  <- n_init*0.4}
  if (c>0.7  & c<=0.8 )       {inflation_f   <- 1.5  ; min.opt <- n_init*0.7}
  if (c>0.8  & c<=0.85)       {inflation_f   <- 2.1    ; min.opt <- n_init*0.8}
  if (c>0.85 & c<=0.9)        {inflation_f   <- 2.8  ; min.opt <- n_init*0.9}


  if (c<=0.7  & n.predictors <6)            {inflation_f    <- 0.9 ; min.opt  <- n_init*0.3}
  if (c>0.7  & c<=0.8  & n.predictors < 8)   {inflation_f   <- 1.8  ; min.opt <- n_init*0.4}
  if (c>0.8  & c<=0.85 & n.predictors < 8)   {inflation_f   <- 2.1    ; min.opt <- n_init *0.5}
  if (c>0.85 & c<=0.9  & n.predictors < 8)   {inflation_f   <- 2.8  ; min.opt <- n_init*0.7}

  max.opt <- inflation_f * n_init

  tol = max(5,ceiling(round(n_init/200)/5) * 5)

  print("Optimisation Starting, ~ 1 min left...")

  #Automatically adjust number of simulations to ensure MCSE is not too high
  A   <- 2*p*(1-p)*stats::qnorm(c)^2
  app <- sqrt(1/(A* max.opt)+2/(max.opt-2) )

  # if (app/sqrt(nsim)>0.0027) nsim = ceiling(app^2/0.0025^2/100)*100

  # s_est_quick <- function(n, nsim=100){
  #
  #   s <-  expected_s_n_binary_quick(n, S = S, mean_eta = mean_eta, variance_eta = variance_eta,  beta = beta, p = p, c = c, n.predictors = n.predictors, nval = nval, nsim = 100,  parallel=parallel)
  #   #(round(s[1]/0.0025)*0.0025-s[2]) - S
  #   s[1] - S
  # }
  #
  # n <- bisection(s_est_quick, min.opt, max.opt, tol = tol, nsim = 100)
  # tol = ceiling(round(n_init/200)/10) * 10
  # n <- ceiling(n/tol)*tol
  #
  # max.opt <- n*1.15
  # min.opt <- n*0.9


  s_est <- function(n, nsim=nsim){

    s <-  expected_s_n_binary(n, S = S, mean_eta = mean_eta, variance_eta = variance_eta,  beta = beta, p = p, c = c, n.predictors = n.predictors, nval = nval, nsim = nsim, parallel=parallel)
    #(round(s[1]/0.0025)*0.0025-s[2]) - S
    s[1] - S
  }

  n   <- bisection(s_est, min.opt, max.opt, tol = tol, nsim = nsim)
  tol <- ceiling(round(n_init/200)/5) * 5
  if (tol==0) tol <- 5
  n   <- ceiling(n/tol)*tol

#   # run <- expected_s(n, p=p, c=c, n.true=n.true, n.noise=n.noise, beta = c(0.5,0.3,0.3,0.15,0.15), nsim=1000, nval=50000, cores=2)

  #print(paste("Required sample size: ", n ))

  size               <- NULL
  size$rvs           <- as.vector(n_init)
  size$sim           <- as.vector(n)
  # size$n_simulations <- nsim
  # size$correct_to_nearest <- as.vector(tol)

  size

}
