
#' Sample size required to develop a risk prediction model for binary outcomes
#'
#' @description
#' This function calculates the sample size required to achieve an expected Calibration Slope (S), given anticipated features of the data and the model
#' (outcome prevalence, C-statistic and number of predictors).
#'
#'It takes approximately one minute to run. Ideally it should be followed by checking also
#' the Mean Absolute Prediction Error that corresponds to the calculated sample size.

#' @param MAPE (numeric) The target expected calibration slope
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


samplesizedev_binary_mape_corr <- function(MAPE, p, c,  n.predictors, beta, nval = 25000, nsim = 1000, parallel = TRUE, cor0=0, cor1=0){

  set.seed(2022)

  mean_var         <- find_mu_sigma(p,c)
  mean_eta         <- mean_var[1]
  variance_eta     <- mean_var[2]

  # Find beta that corresponds to that variance

  if (cor0==0 & cor1 ==0) {

    betan    <- beta * sqrt(mean_var[2]/sum(beta^2))
    sigma   <- diag(1, n.predictors)} else

    {

      betan  <- adjust_multiplier_correlated(c=c, mean = mean_eta, beta = beta, n.predictors = n.predictors, cor0=cor0, cor1=cor1)

      n.noise <- length(beta[beta==0])
      n.true  <- n.predictors-n.noise


      # Specify correlation matrix
      sigma <- matrix(0, nrow = n.predictors,  ncol = n.predictors)
      sigma[1:n.true, 1:n.true] <- cor0
      if (n.noise>0) {
        sigma[(n.true+1):n.predictors, (n.true+1):n.predictors] <- cor1}
      diag(sigma) <- 1

    }

  #r2   <- as.numeric(approximate_R2(c, p, n = 1000000)[2])

  n_init <- exp((-0.508 + 0.259 * log(p) + 0.504 * log(n.predictors) - log(MAPE))/0.544) ;

  min.opt = round(n_init*0.5)
  max.opt = round(n_init*1.5)


  tol = ceiling(round(n_init/200)/5) * 5
  #tol = 20

  print("Optimisation Starting ~ 1 min left...")

  mape_est <- function(n, nsim=nsim){

    mape <-  expected_mape_n_binary_corr(n, MAPE = MAPE, mean_eta = mean_eta, variance_eta = variance_eta, beta=betan,  p = p, c = c, n.predictors = n.predictors, nval = nval, nsim = nsim, parallel = parallel, cor0=cor0)

    MAPE-mape[1]

  }


  n <- bisection_mape(mape_est, MAPE=MAPE, min.opt, max.opt, tol = tol, nsim = nsim)
  tol = ceiling(round(n/200)/5) * 5
  n <- ceiling(n/tol)*tol

  #run <- expected_s(n, p=p, c=c, n.true=n.true, n.noise=n.noise, beta = c(0.5,0.3,0.3,0.15,0.15), nsim=1000, nval=50000, cores=2)

  #print(paste("Required sample size: ", n ))

  size        <- NULL
  size$rvs2   <- as.vector(round(n_init))
  size$actual <- as.vector(round(n))

  size

}
