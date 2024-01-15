
#' Sample Size Required to Develop a Risk Prediction model for Binary Outcomes
#'
#' @description
#' This function calculates the sample size required to achieve an expected Calibration Slope (S), given anticipated features of the data and the model. For a binary outcome these
#' features are the outcome prevalence (phi) and c-statistic (c) and the number of candidate predictor variables (p).
#' Future versions will  accept the anticipated proportion of events at a given time point (phi) and c-index (c)  for a Survival outcome.
#'
#' The required value for the (average) calibration slope should be at least S=0.9. On the other hand, the input
#' values for phi, c, and p are study-dependent.
#'
#'
#' The calculation takes usually takes 1-2 minutes for a binary outcome.
#'
#' The variability in the calibration slope, can be seen as an indicator of model stability
#' for the calculated sample size. This can be evaluated using the function 'expected_cs'. For Binary
#' outcomes, this check also provides the expected Mean Absolute Prediction Error (MAPE) that corresponds
#' to the calculated sample size.

#' @param outcome (character) The type of outcome (''Binary'' ; ''Survival'' to be added in later versions)
#' @param S (numeric) The target expected calibration slope
#' @param phi (numeric) The anticipated outcome prevalence of the binary outcome (or proportion of events for survival outcome)
#' @param c (numeric) The anticipated c-statistic for binary outcome (or c-index for survival outcome)
#' @param p (numeric) The number of candidate predictor variables
#' @param gamma (vector) The relative strength of predictors (0 for noise)
#' @param nsim (numeric) The number of simulations (at least 500, default value 1000 to ensure small simulation error)
#' @param nval (numeric) Size of validation data (at least 10000)
#' @param parallel (logical) parallel processing to speed up computations (default=TRUE)

#'
#' @return n: the required sample size
#' @export
#'
#' @examples
#' # Binary Outcome: Find the sample size required for an average calibration slope of S = 0.9
#' # samplesizedev(outcome="Binary", S = 0.9, phi = 0.2, c = 0.85, p = 10, parallel = FALSE)
#'
#' # Binary Outcome: Prefer parallel computing  that ensures faster run
#' # samplesizedev(outcome="Binary", S = 0.9, phi = 0.2, c = 0.85, p = 10, parallel=TRUE)
#'
#' # Binary Outcome: Check the expected MAPE and Calibration Slope for the selected size
#' # expected_cs (outcome= "Binary", n = 530, phi= 0.2, c = 0.85, p = 10, nsim = 1000)
#'
#' # Survival Outcome: Find the sample size required for an average calibration slope of S = 0.9
#' # samplesizedev(outcome = "Survival", S = 0.9, phi = 0.2, c = 0.85, p = 10,  nsim = 500, parallel = TRUE)
#'
#' # Survival  Outcome: Check the expected MAPE and Calibration Slope for the selected size
#' # expected_cs(outcome = "Survival", n = 390, phi = 0.2, c = 0.85, p = 10, nsim = 500, parallel = TRUE)
#'
#'
#' @seealso
#' expected_cs


samplesizedev <- function(outcome="Binary", S = NULL, MAPE = NULL, phi, c,  p, gamma = rep(1/p, p), nval = 25000, nsim = 1000, parallel = TRUE){

  beta          <- gamma
  n.predictors  <- p
  p             <- phi

  if (outcome=="Binary")   { if (length(MAPE)==0)       n <- samplesizedev_binary_s(S=S, p=p, c=c,  n.predictors = n.predictors, beta=beta, nval = nval, nsim = nsim, parallel = parallel) else
                             if (length(S)==0)          n <- samplesizedev_binary_mape(MAPE=MAPE, p=p, c=c,  n.predictors = n.predictors, beta=beta, nval = nval, nsim = nsim, parallel = parallel)}

  if (outcome=="Survival") nval = 10000
  if (outcome=="Survival") n <- samplesizedev_survival(S=S, p=p, c=c,   n.predictors = n.predictors, beta=beta, nval = nval, nsim = nsim, parallel = parallel)

  n

}
