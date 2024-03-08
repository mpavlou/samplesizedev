
<!-- README.md is generated from README.Rmd. Please edit that file -->

# samplesizedev

<!-- badges: start -->
<!-- badges: end -->

The package samplesizedev performs unbiased sample size calculations
(using simulation) for the development of risk models for binary
outcomes. It requires information on the anticipate values of the:

- outcome prevalence

- c-statistic (AUC)

- number of predictor variables

to calculate the sample size required to achieve an calibration slope
(S) or Mean Absolute Prediction Error (MAPE), on average.

## Installation

The development version of samplesizedev can be installed from
[GitHub](https://github.com/) with:

``` r
# If package 'devtools is not installed, first install it'
# install.packages("devtools")
# require("devtools")

devtools::install_github("mpavlou/samplesizedev")
require(samplesizedev)
```

## Example

This is a basic example which shows how to:

1)  to calculate the expected calibration slope and MAPE for a given
    sample size, and

2)  To calculate the sample size to achieve a target expected
    calibration slope of S=0.9 for a binary outcome

``` r
library(samplesizedev)

# Explore the two main commands:
# ?samplesizedev
# ?expected_cs

# Calculate the expected calibration slope and MAPE
# Sample size=500; Prevalence=0.2; C-statistic=0.8; Number of predictors=10; 


expected_cs(n = 500, phi = 0.2, c = 0.85, p = 10)
```

<img src="man/figures/README-example-1.png" width="100%" />

    #>     N Mean_CS  SD_CS RMSD_CS Pr(CS<0.8) Mean_MAPE SD_MAPE Prev. C-Stat.
    #> 1 500   0.902 0.1002  0.1404       0.16    0.0393  0.0087   0.2    0.85
    #>    # Predictors
    #> 1            10

    # Calculate sample size for target calibration slope
    # Target Calibration slope S=0.9; Prevalence=0.2; c-statistic=0.8; Number of predictors=10; 

    samplesizedev(outcome="Binary", S = 0.9, phi = 0.2, c = 0.85, p= 10)
    #> [1] "Optimisation Starting, ~ 1-2 min left..."
    #> $rvs
    #> [1] 308
    #> 
    #> $sim
    #> [1] 500
