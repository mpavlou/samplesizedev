
#' Calculate the expected calibration slope, MAPE and other performance metrics for a given  sample size
#'
#' @description
#' This function calculates the expected calibration slope and MAPE given key data and model characteristics
#' (outcome prevalence, c-statistic and number of predictors). It takes approximately 15 seconds to run (binary outcome).
#'
#' @param outcome (character) The type of outcome (''Binary'' or ''Survival'')
#' @param n (numeric) The sample size
#' @param phi (numeric) The anticipated outcome prevalence (binary outcome) or proportion of events (survival outcome)
#' @param c (numeric) The anticipated c-statistic (binary outcome) or c-index (survival outcome)
#' @param p (numeric) The number of candidate predictor variables
#' @param nsim (numeric) The number of simulations (default=1000; use at least 500 to ensure small simulation error)
#' @param nval (numeric) Size of validation data (default=25000; use at least 10000)
#' @param parallel (logical) parallel processing to speed up computations (default=TRUE)
#' @param method (character) the fitting method (default="MLE"; currently the only option. Other fitting methods will be added in future versions)
#' @param gamma (numeric) Relative strength of predictors (default=rep(1/p,p); same length as p, must sums up to 1)
#' @param long (logical) Extract results for all simulations instead of just averages

#'
#' @return   a data frame df with the following elements:
#' @return   input sample size (n)
#' @return   expected calibration slope (mean_CS)
#' @return   standard deviation of the CS (sd_CS)
#' @return   the probability of obtaining a mis-calibrated model with calibration slope <0.8 (Pr(CS<0.8))
#' @return   the expected MAPE (MAPE)
#' @return   the standard deviation of the expected MAPE (sd_MAPE)
#' @return   the expected optimism in R square Nagelgerke (optimism_R2_Nag)
#'
#'
#' @export
#'
#' @examples
#' # expected_performance(outcome="Binary", n = 530, phi = 0.2, c = 0.85, p = 10, nsim = 100, parallel = FALSE)
#'
#' # Prefer parallel computing with >2 cores that ensure faster running
#' # expected_performance(n = 530, phi = 0.2, c = 0.85, p = 10, nsim = 100, parallel = TRUE)

#' @seealso
#' samplesizedev

#'
#'
expected_performance <- function(outcome="Binary", n, phi, c,  p,  gamma = rep(1/p, p), nval = 25000, nsim = 1000, parallel = TRUE, method = "MLE", long=FALSE, approx = FALSE, x=NULL, y=NULL){

  beta          <- gamma
  n.predictors  <- p
  p             <- phi

  if (outcome=="Binary" & length(x)==0)   performance <- expected_cs_mape_binary (n=n,  p=p, c=c,  n.predictors = n.predictors, beta=beta, nval = nval, nsim = nsim, parallel = parallel, method = method, long = long, approx = approx)
  if (outcome=="Binary" & length(x)!=0)   performance <- expected_cs_mape_binary_xy (n=n,  p=p, c=c,  n.predictors = n.predictors, beta=beta, nval = nval, nsim = nsim, parallel = parallel, method = method, long = long, approx = approx, x=x, y=y)

  if (outcome=="Survival") performance <- expected_cs_survival (n=n, p=p, c=c,  n.predictors = n.predictors, nval = nval, nsim = nsim, parallel = parallel, method = method)

  performance

}
