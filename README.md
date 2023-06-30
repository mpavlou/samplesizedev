
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

This is a basic example which shows a) how to calculate the expected
shrinkage and MAPE for a given sample size, and b) To calculate the
sample size to achieve a taregt expected shrinkage of S=0.9

``` r
library(samplesizedev)

expected_cs(n = 530, p = 0.2, c = 0.85, n.predictors = 10, nsim = 100, parallel = TRUE)
```

<img src="man/figures/README-example-1.png" width="100%" />

    #>     N Expected CS SD(CS) RMSD(CS) Pr(CS<0.8) Expected MAPE SD(MAPE) Prevalence
    #> 1 530         0.9 0.0927   0.1357       0.15         0.038   0.0088        0.2
    #>   C-Statistic  # Predictors
    #> 1        0.85            10

    samplesizedev(outcome="Binary", S = 0.9, p = 0.2, c = 0.85, n.predictors = 10,  nsim = 1000)
    #> [1] "Optimisation Starting ~ 1 min left..."
    #> $riley
    #> [1] 305
    #> 
    #> $actual
    #> [1] 530

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
