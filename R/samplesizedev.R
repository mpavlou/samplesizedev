
#' Sample Size Required to Develop a Risk Prediction model for Binary and Survival outcomes
#'
#' @description
#' This function calculates the sample size required to achieve an expected Calibration Slope (S), given anticipated features of the data and the model. These
#' features are the outcome prevalence (p) and c-statistic (c) for a Binary outcome  or the proportion of events at a given time point (p) and c-index (c)  for a survival outcome,
#' and the number of candidate predictor variables (n.predictors).
#'
#' The required value for the (average) calibration slope should be at least S=0.9. On the other hand, the input
#' values for p, c, and n.predictors are application-dependent.
#'
#'
#' The calculation takes approximately one minute for binary outcomes and 2-4 minutes for survival outcome (depending on the number of predictors).
#'
#' We suggest that the variability in the calibration
#' slope, which is an indicator of model stability for the calculated sample size, be subsequently checked using the function 'expected_cs'. For binary outcomes, this
#' check also provides the expected Mean Absolute Prediction Error (MAPE) that corresponds to the calculated sample size.

#' @param outcome (character) The type of outcome (''Binary'' or ''Survival'')
#' @param S (numeric) The target expected calibration slope
#' @param p (numeric) The anticipated outcome prevalence (binary outcome) or proportion of events (survival outcome)
#' @param c (numeric) The anticipated c-statistic (binary outcome) or c-Index (survival outcome)
#' @param n.predictors (numeric) The number of candidate predictor variables
#' @param beta (vector) The relative strength of predictors (0 for noise)
#' @param nsim (numeric) The number of simulations (at least 500, default value 1000 to ensure small simulation error)
#' @param nval (numeric) Size of validation data (at least 10000)
#' @param parallel (logical) parallel processing to speed up computations (default=TRUE)

#'
#' @return n: the required sample size
#' @export
#'
#' @examples
#' # Binary Outcome: Find the sample size required for an average calibration slope of S = 0.9
#' # samplesizedev(outcome="Binary", S = 0.9, phi = 0.2, c = 0.85, p = 10,  nsim = 500, parallel = FALSE)
#'
#' # Binary Outcome: Prefer parallel computing with >2 cores that ensure faster running
#' # samplesizedev(outcome="Binary", S = 0.9, phi = 0.2, c = 0.85, p = 10,  nsim = 500, parallel = TRUE)
#'
#' # Binary Outcome: Check the expected MAPE and Calibration Slope for the selected size
#' # expected_cs (outcome= "Binary", n = 530, phi= 0.2, c = 0.85, p = 10, nsim = 500, parallel = TRUE)
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

  if (outcome=="Survival") n <- samplesizedev_survival(S=S, p=p, c=c,   n.predictors = n.predictors, beta=beta, nval = nval, nsim = nsim, parallel = parallel)

  n

}
