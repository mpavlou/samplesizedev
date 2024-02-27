
<!-- README.md is generated from README.Rmd. Please edit that file -->

# samplesizedev

<!-- badges: start -->
<!-- badges: end -->

The package samplesizedev performs unbiased sample size calculations for
the development of risk models.

## Installation

You can install the development version of samplesizedev from
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

1)  to calculate the expected shrinkage and MAPE for a given sample
    size, and

2)  To calculate the sample size to achieve a target expected shrinkage
    of S=0.9 for a binary outcome

``` r
library(samplesizedev)

# Sample size=500; Prevalence=0.2; cC-statistic=0.8; Number of predictors=10; 
# Calculate expected calibration slope and MAPE

expected_cs(n = 500, phi = 0.2, c = 0.85, p = 10, nsim = 1000)
```

<img src="man/figures/README-example-1.png" width="100%" />

    #>     N Mean_CS  SD_CS Pr(CS<0.8) Mean_MAPE SD_MAPE Prev. C-Stat.  # Predictors
    #> 1 500   0.902 0.1002       0.16    0.0393  0.0087   0.2    0.85            10


    # Target Calibration slope=0.9; Prevalence=0.2; c-statistic=0.8; Number of predictors=10; 
    # Calculate sample size

    samplesizedev(outcome="Binary", S = 0.9, phi = 0.2, c = 0.85, p= 10,  nsim = 1000)
    #> [1] "Optimisation Starting, ~ 1-2 min left..."
    #> $rvs
    #> [1] 308
    #> 
    #> $sim
    #> [1] 500
