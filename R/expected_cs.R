
#' Calculate the expected calibration slope (and mean absolute prediction error (MAPE) for binary data) for a given  sample size
#'
#' @description
#' This function calculates the expected calibration slope and MAPE given key data and model characteristics
#' (outcome prevalence/, C-statistic and number of predictors). It takes approximately 15 seconds to run (binary outcome).
#'
#' @param outcome (character) The type of outcome (''Binary'' or ''Survival'')
#' @param n (numeric) The sample size
#' @param p (numeric) The anticipated outcome prevalence (binary outcome) or proportion of events (survival outcome)
#' @param c (numeric) The anticipated C-statistic (binary outcome) or C-Index (survival outcome)
#' @param n.predictors (numeric) The number of candidate predictor variables
#' @param nsim (numeric) The number of simulations (at least 500, default value 1000 to ensure small simulation error)
#' @param nval (numeric) Size of validation data (at least 10000)
#' @param parallel (logical) parallel processing to speed up computations (default=TRUE)
#' @param method (character) the fitting method. "MLE" is the default and currently only option, but others will be added in future versions
#' @param beta (numeric) the relative strength of predcitors (same lenght as n.predictors)


#'
#' @return df: the expected calibration slope (and MAPE for binary outcomes)
#' @export
#'
#' @examples
#' # expected_cs(outcome="Binary", n = 530, p = 0.2, c = 0.85, n.predictors = 10, nsim = 100, parallel = FALSE)
#'
#' # Prefer parallel computing with >2 cores that ensure faster running
#' # expected_cs_mape_binary(n = 530, p = 0.2, c = 0.85, n.predictors = 10, nsim = 100, parallel = TRUE)

#' @seealso
#' samplesizedev

#'
#'
expected_cs <- function(outcome="Binary", n, p, c,  n.predictors, nval = 25000, nsim = 1000, parallel = TRUE, method = "MLE"){

  if (outcome=="Binary")   performance <- expected_cs_mape_binary (n=n,  p=p, c=c,  n.predictors = n.predictors, nval = nval, nsim = nsim, parallel = parallel, method = method)

  if (outcome=="Survival") performance <- expected_cs_survival (n=n, p=p, c=c,  n.predictors = n.predictors, nval = nval, nsim = nsim, parallel = parallel, method = method)

  performance

}
