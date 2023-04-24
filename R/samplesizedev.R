
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
#' @export
#'
#' @examples
#' # Find the sample size
#'   samplesizedev(S = 0.9, p = 0.2, c = 0.85, n.predictors = 10,  nsim = 500, parallel = FALSE)
#'
#' # Prefer parallel computing with >2 cores that ensure faster running
#' # samplesizedev(S = 0.9, p = 0.2, c = 0.85, n.predictors = 10,  nsim = 1000, parallel = TRUE)
#'
#' # Check the expected MAPE and Calibration Slope for the selected size
#' # expected_cs_mape(n = 530, p = 0.2, c = 0.85, n.predictors = 10, nsim = 1000, parallel = TRUE)
#'
#'
#'
#' @seealso
#' expected_cs_mape


samplesizedev <- function(S, p, c,  n.predictors, nval = 25000, nsim = 1000, parallel = TRUE){

  set.seed(1)

  mean_var_eta     <- find_mu_sigma(p, c)
  mean_eta         <- mean_var_eta[1]
  variance_eta     <- mean_var_eta[2]

  r2   <- as.numeric(approximate_R2(c, p, n = 100000)[2])

  n_init <- round((n.predictors)/ ((S-1)*log(1-r2/S)))

  min.opt                              <- n_init*0.8
  if (c<=0.7)            inflation_f   <- 1.3
  if (c>0.7  & c<=0.8)   inflation_f   <- 1.5
  if (c>0.8  & c<=0.85)  inflation_f   <- 2
  if (c>0.85 & c<=0.9)   inflation_f   <- 2.8
  max.opt                              <- inflation_f*n_init

  tol = round(n_init/100/10)*10
  #tol = 20

  print("Optimisation Starting ~ 1 min left...")
  s_est <- function(n, nsim=nsim){

    s <-  expected_s_n(round(n), S = S, mean_eta = mean_eta, variance_eta = variance_eta,  p = p, c = c, n.predictors = n.predictors, nval = nval, nsim = nsim, parallel=parallel)[1]
    round(s/0.0025)*0.0025-S
    #(s - S)^2
    #abs(s-S)
  }

  n <- bisection(s_est, min.opt, max.opt, tol = tol, nsim = nsim)
  n <- ceiling(n/10)*10

  #run <- expected_s(n, p=p, c=c, n.true=n.true, n.noise=n.noise, beta = c(0.5,0.3,0.3,0.15,0.15), nsim=1000, nval=50000, cores=2)

  print(paste("Required sample size: ", n ))

}
