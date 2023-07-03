
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
# install.packages("devtools")
devtools::install_github("mpavlou/samplesizedev")
```

## Example

This is a basic example which shows how to a) to calculate the expected
shrinkage and MAPE for a given sample size, and b) To calculate the
sample size to achieve a target expected shrinkage of S=0.9 for a binary
outcome

``` r
library(samplesizedev)

expected_cs(n = 530, p = 0.2, c = 0.85, n.predictors = 10, nsim = 1000, parallel = TRUE)
```

<img src="man/figures/README-example-1.png" width="100%" />

    #>     N Mean_CS  SD_CS RMSD_CS Pr(CS<0.8) Mean_MAPE SD_MAPE Prev. C-Stat.
    #> 1 530     0.9 0.0979   0.139       0.14    0.0384  0.0086   0.2    0.85
    #>    # Predictors
    #> 1            10

    samplesizedev(outcome="Binary", S = 0.9, p = 0.2, c = 0.85, n.predictors = 10,  nsim = 1000)
    #> [1] "Optimisation Starting ~ 1 min left..."
    #> $riley
    #> [1] 305
    #> 
    #> $actual
    #> [1] 530
