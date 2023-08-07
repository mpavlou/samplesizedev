
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


samplesizedev_binary_s_corr <- function(S, p, c, n.predictors, beta=rep(1/n.predictors, n.predictors), nval = 25000, nsim = 1000, parallel = TRUE, cor0=0.1, cor1=0.05){

  set.seed(2022)

  mean_var     <- find_mu_sigma(p, c, tol = 0.00001)
  mean_eta         <- mean_var[1]
  variance_eta     <- mean_var[2]

  # Find beta that corresponds to that variance

  if (cor0==0 & cor1 ==0) {

    betan    <- beta * sqrt(mean_var[2]/sum(beta^2))
    sigma   <- diag(1, n.predictors)} else

    {

      betan  <- adjust_multiplier_correlated(c=c, mean =  mean_eta, beta = beta, n.predictors = n.predictors, cor0=cor0, cor1=cor1)

      n.noise <- length(beta[beta==0])
      n.true  <- n.predictors-n.noise


      # Specify correlation matrix
      sigma <- matrix(0, nrow = n.predictors,  ncol = n.predictors)
      sigma[1:n.true, 1:n.true] <- cor0
      if (n.noise>0) {
        sigma[(n.true+1):n.predictors, (n.true+1):n.predictors] <- cor1}
      diag(sigma) <- 1

    }


  # xval    <- mvtnorm::rmvnorm(1000000, rep(0, n.predictors), sigma = sigma)
  # yval    <- stats::rbinom(1000000, 1, invlogit(mean_eta + xval%*%betan))
  #
  # lp <- mean_eta + xval%*%betan
  #
  # cs           <- glm(yval~lp, family="binomial")
  # cs.ml        <- coef(cs)[2]
  # r2.ml        <- PseudoR2(cs,"CoxSnell")
  # r2.ml
  #
  # n_init <- round((n.predictors)/ ((S-1)*log(1-r2.ml/S)))
  # n_init

  r2   <- as.numeric(approximate_R2(c, p, n = 200000)[2])

  n_init <- round((n.predictors)/ ((S-1)*log(1-r2/S)))



  if (c<=0.7)            {inflation_f   <- 1.3 ; min.opt  <- n_init*0.7}
  if (c>0.7  & c<=0.8)   {inflation_f   <- 1.6  ; min.opt <- n_init}
  if (c>0.8  & c<=0.85)  {inflation_f   <- 2.1    ; min.opt <- n_init}
  if (c>0.85 & c<=0.9)   {inflation_f   <- 2.8  ; min.opt <- n_init}

  max.opt  <- inflation_f*n_init

  tol = ceiling(round(n_init/200)/5) * 5

  print("Optimisation Starting ~ 1 min left...")
  s_est <- function(n, nsim=nsim){

    s <-  expected_s_n_binary_corr(n, S = S, mean_eta = mean_eta, variance_eta = variance_eta,  beta = betan, p = p, c = c, n.predictors = n.predictors, nval = nval, nsim = nsim, parallel=parallel, cor0=cor0, cor1=cor1)
    s[1] - S
  }

  n <- bisection(s_est, min.opt, max.opt, tol = tol, nsim = nsim)

  tol <- ceiling(round(n/200)/5) * 5
  n   <- ceiling(n/tol)*tol

  size        <- NULL
  size$rvs1   <- as.vector(round(n_init))
  size$actual <- as.vector(round(n))

  size

}

# Examples

# system.time(n1 <- samplesizedev_binary_s_corr(S=0.9, p = 0.174, c = 0.89, n.predictors = 24, nsim = 2000, parallel = TRUE, nval = 25000))
#
# n1
#
# expected_cs_mape_binary_corr(n =  n1$actual, p = 0.174, c = 0.89, n.predictors = 24, nsim = 2000, parallel = TRUE, nval=25000)
#
#
# n1$rvs1*p/12
#
# n1$actual*p/12
