
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


samplesizedev_binary_mape <- function(MAPE, p, c,  n.predictors, beta, nval = 25000, nsim = 1000, parallel = TRUE){

  if (p>0.5) p  <- 1-p

  set.seed(1)

  mean_var_eta     <- find_mu_sigma(p, c)
  mean_eta         <- mean_var_eta[1]
  variance_eta     <- mean_var_eta[2]

  #r2   <- as.numeric(approximate_R2(c, p, n = 1000000)[2])

  n_init <- exp((-0.508 + 0.259 * log(p) + 0.504 * log(n.predictors) - log(MAPE))/0.544) ;


  min.opt = round(n_init*0.5)
  max.opt = round(n_init*1.5)


  tol = ceiling(round(n_init/200)/10) * 10
  #tol = 20

  print("Optimisation Starting ~ 1 min left...")

  mape_est <- function(n, nsim=nsim){

    mape <-  expected_mape_n_binary(n, MAPE = MAPE, mean_eta = mean_eta, variance_eta = variance_eta,  p = p, c = c, beta = beta, n.predictors = n.predictors, nval = nval, nsim = nsim, parallel = parallel)

    MAPE-round(mape[1]/0.00025)*0.00025
    MAPE-round(mape[1]/0.005)*0.005
    mape[1] - MAPE/phi

  }


  n <- bisection_mape(mape_est, min.opt, max.opt, tol = tol, nsim = nsim)
  tol = ceiling(round(n_init/200)/5) * 5
  n <- ceiling(n/tol)*tol

  #run <- expected_s(n, p=p, c=c, n.true=n.true, n.noise=n.noise, beta = c(0.5,0.3,0.3,0.15,0.15), nsim=1000, nval=50000, cores=2)

  #print(paste("Required sample size: ", n ))

  size        <- NULL
  size$rvs2   <- as.vector(round(n_init))
  size$actual <- as.vector(n)

  size

}
