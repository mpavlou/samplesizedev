#' Calculate the expected calibration slope, MAPE and other performance metrics for a given  sample size
#'
#' @description
#' This function calculates the sampling distribution of calibration slope and MAPE and other metrics given key data and model characteristics
#' (outcome prevalence, c-statistic and number of predictors). It takes approximately 15 seconds to run (binary outcome).
#'
#' @param outcome (character) The type of outcome (''Binary'' or ''Survival'')
#' @param n (numeric) The sample size
#' @param phi (numeric) The anticipated outcome prevalence (binary outcome) or proportion of events (survival outcome)
#' @param c (numeric) The anticipated c-statistic (binary outcome) or c-index (survival outcome)
#' @param p (numeric) The number of candidate predictor variables
#' @param nsim (numeric) The number of simulations (default=1000; use at least 500 to ensure small Monte Carlo simulation error)
#' @param nval (numeric) Size of validation data (default=25000; use at least 10000)
#' @param parallel (logical) parallel processing to speed up computations (default=TRUE)
#' @param method (character) the fitting method (default="MLE"; currently the only option. Other fitting methods will be added in future versions)
#' @param gamma (numeric) Relative strength of predictors (default=rep(1/p,p); same length as p, must sums up to 1)
#' @param long (logical) Extract results for all simulations instead of just averages
#' @param threshold (numeric) A threshold to calculate sensitivity (default=outcome prevalence)
#' @param individual_predicted_probability (numeric) An individual predicted probability of inderest (default at median predicted probability) to obtain the uncertainty due to finite sample size


#' @return   Data frame df with the following elements:
#' @return   Input sample size (n)
#' @return   Specified outcome prevalence (phi)
#' @return   Specified C-statistic / C-index
#' @return   Specified Number of predictors
#' @return   Expected calibration slope (mean_CS)
#' @return   Standard deviation of the CS (sd_CS)
#' @return   Probability of obtaining a model with calibration slope between 0.85 and 1.15
#' @return   Expected MAPE (MAPE)
#' @return   Standard deviation of the expected MAPE (sd_MAPE)
#' @return   Expected optimism in R square Nagelgerke (optimism_R2_Nag)
#' @return   Expected AUC
#' @return   Standard deviation of the expected AUC
#' @return   Standard deviation of the average predicted risk
#' @return   Expected Brier score
#' @return   Sensitivity at the specified threshold (Sensitivity)
#' @return   Net-Benefit at the specified threshold
#' @return   Expected Individual predicted risk at the specified percentile
#' @return   Standard deviation of the Individual predicted probability (IPP)
#'
#'
#' @import foreach
#' @importFrom foreach %dopar%
#' @import doParallel
#' @import doRNG
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
expected_performance <- function(outcome="Binary", n, phi, c,  p,  gamma = rep(1/p, p), nval = 25000, nsim = 1000, parallel = TRUE, method = "MLE", long=FALSE, approx = FALSE, x=NULL, y=NULL, threshold=phi, individual_predicted_probability = NULL){

  beta          <- gamma
  n.predictors  <- p
  p             <- phi

  if (outcome=="Binary" & length(x)==0)   performance <- expected_cs_mape_binary (n=n,  p=p, c=c,  n.predictors = n.predictors, beta=beta, nval = nval, nsim = nsim, parallel = parallel, method = method, long = long, approx = approx, threshold = threshold, individual_predicted_probability = individual_predicted_probability)
  if (outcome=="Binary" & length(x)!=0)   performance <- expected_cs_mape_binary_xy (n=n,  p=p, c=c,  n.predictors = n.predictors, beta=beta, nval = nval, nsim = nsim, parallel = parallel, method = method, long = long, approx = approx, x=x, y=y, threshold = threshold, individual_predicted_probability = individual_predicted_probability)

  if (outcome=="Survival") performance <- expected_cs_survival (n=n, p=p, c=c,  n.predictors = n.predictors, nval = nval, nsim = nsim, parallel = parallel, method = method)

  performance

}
