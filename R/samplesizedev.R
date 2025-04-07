
#' Sample Size Required to Develop a Risk Prediction model for Binary Outcomes
#'
#' @description
#' This function calculates the sample size required to achieve an expected Calibration Slope (S) or an expected MAPE (MAPE) given anticipated features of the data and the model. For a binary outcome these
#' features are the outcome prevalence (phi) and c-statistic (c) and the number of candidate predictor variables (p).
#' Future versions will  accept the anticipated proportion of events at a given time point (phi) and c-index (c)  for a time to event outcome.
#'
#' The user first needs to specify the study-depended input values for the outcome prevalence the c-statistic and the number of predictors.
#'
#' The user also needs to input a target value for *either* the (expected) calibration slope (S) or the (expected) mean absolute prediction error (MAPE).
#' The value of S should be at least 0.9, to ensure that the degree of overfitting will be small.  On the other hand, a suitable target value for the (average) MAPE is study-depended and linked to outcome
#' prevalence; target MAPE between (prevalence/10 and prevalence/5) can be a reasonable choice in many cases.
#'
#' The calculation usually takes  a minute or less for a binary outcome.
#'
#' The variability in the calibration slope, can be seen as an indicator of model stability
#' for the calculated sample size. This can be evaluated using the function 'expected_cs'. For Binary
#' outcomes, this check also provides the expected Mean Absolute Prediction Error (MAPE) that corresponds
#' to the calculated sample size.

#' @param outcome (character) The type of outcome (''Binary'' ; ''Survival'' to be added in later versions)
#' @param S (numeric) The target expected calibration slope
#' @param MAPE (numeric) The target expected mean absolute prediction error
#' @param phi (numeric) The anticipated outcome prevalence of the binary outcome (or proportion of events for survival outcome)
#' @param c (numeric) The anticipated c-statistic for binary outcome (or c-index for survival outcome)
#' @param p (numeric) The number of candidate predictor variables
#' @param gamma (numeric) The Relative strength of predictors (default=rep(1/p,p); same length as p, must sum up to 1)
#' @param nsim (numeric) The number of simulations (default=1000; use at least 500 to ensure small simulation error)
#' @param nval (numeric) Size of validation data (default=25000; use at least 10000)
#' @param parallel (logical) parallel processing to speed up computations (default=TRUE)

#'
#' @return n$sim: the required sample size by simulation
#' @return n$rvs: the required sample size by the RvS formulae (only for comparison)
#' @return n$r2_cs: Cox-Snell R square for the input values of phi and p

#' @export
#'
#' @examples
#'
#' # Binary Outcome
#' # Size for target S=0.9
#' # samplesizedev(outcome = "Binary", S = 0.9, phi = 0.2, c = 0.85, p = 10, parallel=TRUE)
#'
#'#' Size for target MAPE=0.04
#' # samplesizedev(outcome = "Binary", MAPE = 0.04, phi = 0.2, c = 0.85, p = 10)
#'
#' # Binary Outcome: Check the expected MAPE and Calibration Slope for a given sample size
#' # expected_performance(outcome = "Binary", n = 530, phi= 0.2, c = 0.85, p = 10)
#'
#' @seealso
#' expected_cs


samplesizedev <- function(outcome="Binary", S = NULL, MAPE = NULL, S1 = NULL, S2 = NULL, P_S1S2 = NULL, phi, c,  p, gamma = rep(1/p, p), nval = 25000, nsim = 1000, parallel = TRUE){

  beta          <- gamma
  n.predictors  <- p
  p             <- phi

<<<<<<< HEAD
  if (n.predictors<=6) nsim = 3000
=======
  if (n.predictors<=6) nsim=3000
>>>>>>> b63b6f2b682f01f907bc7b6a067f2aac0032fca7

  if (outcome=="Binary")   { if (length(MAPE)==0 & length(S1)==0)    n <- samplesizedev_binary_s(S=S, p=p, c=c,  n.predictors = n.predictors, beta=beta, nval = nval, nsim = nsim, parallel = parallel) else
                             if (length(S)==0    & length(S1)==0)    n <- samplesizedev_binary_mape(MAPE=MAPE, p=p, c=c,  n.predictors = n.predictors, beta = beta, nval = nval, nsim = nsim, parallel = parallel) else
                             if (length(S)==0    & length(MAPE)==0)  n <- samplesizedev_binary_prob_s(S1=S1, S2=S2, P_S1S2=P_S1S2, p=p, c=c,  n.predictors = n.predictors, beta = beta, nval = nval, nsim = nsim, parallel = parallel) }

  if (outcome=="Survival") n <- samplesizedev_survival(S=S, p=p, c=c, n.predictors = n.predictors, beta=beta, nval = nval, nsim = nsim, parallel = parallel)

  n

}
