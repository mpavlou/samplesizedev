
#' Sample size required to develop a risk prediction model for binary outcomes
#'
#' @description
#' This function calculates the sample size required to achieve an expected Calibration Slope (S), given anticipated features of the data and the model
#' (outcome prevalence, C-statistic and number of predictors).
#'
#'It takes approximately one minute to run. Ideally it should be followed by checking also
#' the Mean Absolute Prediction Error that corresponds to the calculated sample size.

#' @param l_s (numeric) The lower bound for the calibration slope
#' @param u_s (numeric) The upper bound for the expected calibration slope
#' @param PAP_s (numeric) The probability of acceptable performance in terms of calibration
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
#' expected_performance


samplesizedev_binary_prob_s <- function(l_s, u_s, PAP_s, p, c,   n.predictors, beta = rep(1/n.predictors, n.predictors), nval = 25000, nsim = 1000, parallel = TRUE, plot = TRUE, quick){

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


  a      <- find_n_prap_s(c, p, mean_eta, variance_eta, n.predictors, r2, l_s=l_s, u_s = u_s,
                      PAP_s = PAP_s,min.opt = 0.1, max.opt = 0.99)
  n_init <- a[1]

  if (c <  0.75 )               {inflation_f   <- 1.02  ; min.opt <- n_init *0.95}
  if (c >= 0.75 & c <  0.8  )   {inflation_f   <- 1.1   ; min.opt <- n_init *1}
  if (c >= 0.8  & c <= 0.85 )   {inflation_f   <- 1.2   ; min.opt <- n_init*1.1 }
  if (c >  0.85 & c <= 0.9  )   {inflation_f   <- 1.3   ; min.opt <- n_init*1.2 }

  max.opt <- round(inflation_f * n_init)
  min.opt <- round(min.opt)

  tol = max(5,ceiling(round(n_init/200)/5) * 5)

  if (plot==TRUE & quick==FALSE) print("Optimisation Started: check progress on the appearing plots...") else
    print("Optimisation Started: seconds remaining...")

  #Automatically adjust number of simulations to ensure MCSE is not too high
  A   <- 2*p*(1-p)*stats::qnorm(c)^2
  app <- sqrt(1/(A* max.opt)+2/(max.opt-2) )

  prob_s_est <- function(n, nsim=nsim){

    prob_s <-  expected_prob_s_n_binary(n, l_s = l_s, u_s = u_s, PAP_s = PAP_s, mean_eta = mean_eta,
                                        variance_eta = variance_eta,  beta = beta, p = p, c = c,
                                        n.predictors = n.predictors, nval = nval, nsim = nsim,
                                        parallel=parallel, plot = plot)
    #(round(s[1]/0.0025)*0.0025-s[2]) - S
    prob_s[1] - PAP_s
  }

  n_analytical  <- a[1]

  inflation_f   <- 1.00  ; nc <-   n_analytical*inflation_f

  if (c>0.64  & c<=0.65 )       {inflation_f   <- 0.92  ; nc <-   n_analytical*inflation_f }
  if (c>0.65  & c<=0.66 )       {inflation_f   <- 0.93  ; nc <-   n_analytical*inflation_f }
  if (c>0.66  & c<=0.67 )       {inflation_f   <- 0.95  ; nc <-   n_analytical*inflation_f }
  if (c>0.67  & c<=0.68 )       {inflation_f   <- 0.96  ; nc <-   n_analytical*inflation_f }
  if (c>0.68  & c<=0.69 )       {inflation_f   <- 0.96  ; nc <-   n_analytical*inflation_f }
  if (c>0.69  & c<=0.70 )       {inflation_f   <- 0.96  ; nc <-   n_analytical*inflation_f }
  if (c>0.70  & c<=0.71 )       {inflation_f   <- 0.97  ; nc <-   n_analytical*inflation_f }
  if (c>0.71  & c<=0.72 )       {inflation_f   <- 0.97  ; nc <-   n_analytical*inflation_f }
  if (c>0.72  & c<=0.73 )       {inflation_f   <- 0.99  ; nc <-   n_analytical*inflation_f }
  if (c>0.73  & c<=0.80 )       {inflation_f   <- 1.00  ; nc <-   n_analytical*inflation_f }
  if (c>0.8   & c<=0.81 )       {inflation_f   <- 1.03  ; nc <-   n_analytical*inflation_f }
  if (c>0.81  & c<=0.82 )       {inflation_f   <- 1.14  ; nc <-   n_analytical*inflation_f }
  if (c>0.82  & c<=0.83 )       {inflation_f   <- 1.15  ; nc <-   n_analytical*inflation_f }
  if (c>0.83  & c<=0.84 )       {inflation_f   <- 1.16  ; nc <-   n_analytical*inflation_f }
  if (c>0.84  & c<=0.85 )       {inflation_f   <- 1.08  ; nc <-   n_analytical*inflation_f }
  if (c>0.85  & c<=0.86 )       {inflation_f   <- 1.10  ; nc <-   n_analytical*inflation_f }
  if (c>0.86  & c<=0.87 )       {inflation_f   <- 1.13  ; nc <-   n_analytical*inflation_f }
  if (c>0.87  & c<=0.88 )       {inflation_f   <- 1.19  ; nc <-   n_analytical*inflation_f }
  if (c>0.88  & c<=0.89 )       {inflation_f   <- 1.22  ; nc <-   n_analytical*inflation_f }
  if (c>0.89  & c<=0.90 )       {inflation_f   <- 1.24  ; nc <-   n_analytical*inflation_f }

  n_analytical_corrected <- nc

  if (quick==FALSE) n_sim   <- bisection_prob_s(prob_s_est, min.opt, max.opt, tol = tol, nsim = nsim)

  if (quick==TRUE) {
    size                        <- NULL
    size$analytical_corrected   <- as.vector(round(n_analytical_corrected))
    size$analytical             <- as.vector(round(n_analytical))
  }

  if (quick==FALSE) {
    size                        <- NULL
    size$sim                    <- as.vector(round(n_sim))
    size$analytical_corrected   <- as.vector(round(n_analytical_corrected))
    size$analytical             <- as.vector(round(n_analytical))

  }

  # size$n_simulations <- nsim
  # size$correct_to_nearest <- as.vector(tol)

  size

}
