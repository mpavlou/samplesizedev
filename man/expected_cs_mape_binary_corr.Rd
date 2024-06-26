% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expected_cs_mape_binary_corr.R
\name{expected_cs_mape_binary_corr}
\alias{expected_cs_mape_binary_corr}
\title{Calculate the expected calibration slope and mean absolute prediction error (MAPE) for a given  sample size}
\usage{
expected_cs_mape_binary_corr(
  n,
  p,
  c,
  beta = rep(1/n.predictors, n.predictors),
  n.predictors,
  nsim = 1000,
  nval = 25000,
  cor0 = 0,
  cor1 = 0,
  method = "MLE",
  parallel = TRUE
)
}
\arguments{
\item{n}{(numeric) The sample size}

\item{p}{(numeric) The anticipated outcome prevalence}

\item{c}{(numeric) The C-statistic}

\item{beta}{(numeric) Strength of predictors (same length as n.predictors)}

\item{n.predictors}{(numeric) The number of candidate predictor variables}

\item{nsim}{(numeric) The number of simulations (at least 500, default value 1000)}

\item{nval}{(numeric) Size of validation data}

\item{cor0}{(numeric) correlation between true predictors}

\item{cor1}{(numeric) Correlation between noise predictors}

\item{method}{(character) the fitting method. "MLE" is the default and currently only option, but others will be added in future versions}

\item{parallel}{(numeric) relative strength of predictor variables (same length as n_predictors)}
}
\value{
df: the expected calibration slope and mape
}
\description{
This function calculates the expected calibration slope and MAPE given key data and model characteristics
(outcome prevalence, C-statistic and number of predictors). It takes approximately 15 seconds to test_allrun.
}
\examples{
# expected_cs_mape_binary(n = 530, p = 0.2, c = 0.85, n.predictors = 10, nsim = 100, parallel = FALSE)

# Prefer parallel computing with >2 cores that ensure faster running
# expected_cs_mape_binary(n = 530, p = 0.2, c = 0.85, n.predictors = 10, nsim = 100, parallel = TRUE)
}
\seealso{
samplesizedev_binary samplesizedev_survival
}
