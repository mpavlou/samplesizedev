
<!-- README.md is generated from README.Rmd. Please edit that file -->

# samplesizedev : Simulation-based sample size calculations for the development of risk models 

<!-- badges: start -->
<!-- badges: end -->


R package to calculate the sample size for the development of risk models for binary outcomes.

Related papers: 

**"An evaluation of sample size requirements for developing risk prediction models with binary outcomes"**
published in BMC Medical Research Methodology https://doi.org/10.1186/s12874-024-02268-5

$\textcolor{#f00}{\large  \textbf{NEW:}}$ **"Sample Size Calculations for the Development of Risk Prediction Models that Account for Performance Variability"**
https://doi.org/10.48550/arXiv.2509.14028

### Why do we need 'samplesizedev'?

Riley et al. (2019) proposed 3 formulae, based on 3 distinct criteria, for calculating the sample size for the development of risk models.
- Criterion 1: control overfitting (target: calibration slope=0.9),
- Criterion 2: control optimism in R2 Nageleherke (target: opt=0.05)
- Criterion 3:  precision in the mean predicted risk (target: precision = 0.05).

The formula which aims to control model overfitting ('calibration' formula - C1) most often gives that highest sample size and our article we focused primarily around this formula. While the calibration formula performed well for models with C-statistic/C-index<0.8, we found that that it substantially ***underestimated*** the sample size when the predictive strength of the model was higher. The sample sizes often needed to be increased by 50% or even doubled to meet the calibration targets.

Hence, we developed the **new package 'samplesizedev'** which performs ***unbiased sample size calculations*** regardless of model strength. Our software uses simulation in the background so calculations can take a bit to run (from 30s to some minutes depending on scenatio and computational power). Currently it can be used for the development of risk models for binary outcomes; functionality for ***time to event outcomes***  is under development. 

### $\textcolor{#f00}{\large  \textbf{UPDATE}}$

The package has now been updated and can perform sample size calculations than also ***control the*** $\textcolor{#f00}{variability}$ ***in the calibration slope, instead of *just* the expected value***. 

This is very important because, as shown in the accompanied paper, the variability in performance can be very high when the number of predictors is small. Therefore, while one may think that a reduced model should be preferred to avoid model overfitting, this can be misleading. That's because even though performance is controlled on average,  e.g., E(S)=0.9, the variability can still be very high. In this context, E(S)=0.9 is interpreted to mean that if one were to collect many datasets of the recommended size, and  validate them on large external dataset, then the calibration slope would be on average around 0.9. However, if the variability is very high, the probability of actually obtaining an individual dataset with calibration slope close to 0.9 might be unacceptably low (see examples in the paper and below). 

Hence, in our more recent work we $\textcolor{#f00}{\text{develop sample size calculations where we aim to control the probability of acceptable performance, rather}}$ $\textcolor{#f00}{ \text{than just performance on average}}$. In a simulation-based framework this approach can be easily implemented for any performance metric and a suitably defined range of acceptable performance. Here and in the paper above, we focused on the calibration slope. In addition to the simulation-based approach, we also derived an approximate analytical calculation that avoids and is very quick. 


### How does 'samplesizedev' work?

The software requires information on the anticipated values of the:
- outcome prevalence
- c-statistic (AUC)
- number of predictor variables

Based on the characteristics above it can perform actions based on two core functions:

1. **Calculate the required sample size** to achieve a target expected calibration slope or Mean Absolute Prediction Error (MAPE) (function **'samplesizedev'**)

2. $\textcolor{#f00}{\textbf{NEW:}}$ **Calculate the required sample size** to achieve a $\textcolor{#f00}{\text{a high  probability  of  a   model   with   acceptable   calibration}}$ (function **'samplesizedev'**)

3. **Calculate the expected model performance** at a given sample size (function **'expected_performance'**)


## Installation

The development version of samplesizedev can be installed from
[GitHub](https://github.com/) with:

``` r
# If package 'devtools is not installed, first install it'
 install.packages("devtools")
 require("devtools")

devtools::install_github("mpavlou/samplesizedev")
require(samplesizedev)
```

Please get in touch (m.pavlou@ucl.ac.uk) for any bugs you spot and/or for suggestions for improvement. 

## Examples

This is a basic example which shows how to calculate:

1.  the **sample size** to achieve a target expected calibration slope (e.g. target expected calibration slope **E(S)=0.9**)  or 

2. $\textcolor{#f00}{\text{NEW:}}$ the **sample size**  to achieve a target probability of acceptable performance in terms of calibration (e.g. Probability of calibration slope $\in (0.85,1.15)$, **PrAP(S)=0.8**}

3.  the **expected calibration slope, MAPE and other performance metrics** for a given sample size

``` r
library(samplesizedev)

# Explore the two main commands:
# ?samplesizedev
# ?expected_performance
```

#### Calculation of sample size for given model characteristics, aiming at $\textcolor{#f00}{ \text{ expected Calibration slope,  E(S)=0.9}}$ 

``` r
# Calculate sample size for target calibration slope
# Performance target: E(S)=0.9; Prevalence=0.2; c-statistic=0.85; Number of predictors=10;
# Calculation takes about a minute 

samplesizedev(outcome = "Binary", S = 0.9, phi = 0.2, c = 0.85, p = 10)
#> [1] "Optimisation Started: check progress on the appearing plots..."
#> $rvs
#> [1] 308
#> 
#> $sim
#> [1] 500

# $sim is the sample size calculated by simulation
# $rvs is the sample size calculated using the approach of Riley et al. (2019) (RvS formula Criterion 1 - overfitting)
```

The sample size calculated using simulation is n$sim=500 which corresponds to CS=0.9. In comparison, 
the sample size using previously proposed formulae is n$rvs=308. According to the findings in our paper
the RvS overfitting formula  underestimates the sample size for high C-statistic. Thus, the expected calibration slope will
be in fact lower than we had aim for this size.  We can verify this using the second command of our package, 'expected_performance'.


#### Calculation of sample size for given model characteristics, aiming for Probability of acceptable calibrtion PrAP(S)=0.8 

``` r
# Calculate the sample size Size for Probability of Acceptable Performance (PAP=0.8),
# where Acceptable Performance is defined $S\in (0.85, 1.15)$
# Performance target: PrAP(S)=0.8; Prevalence=0.2; c-statistic=0.85; Number of predictors=10;
samplesizedev(outcome="Binary", l_s= 0.85, u_s = 1.15, PAP_s = 0.8, phi = 0.2, c = 0.85, p = 10)

$sim
[1] 699

# $sim is the sample size calculated by simulation to ensure that PrAP(S)=0.8
```

The sample size calculated using simulation targetting at E(S)=0.9 is 500, while the sample size to ensure that PrAP(S)=0.8 is 699.


#### Calculation of expected model performance (CS, C-statitistc MAPE etc and their variability) for a given sample size and model characteristics

``` r
# Calculate the expected calibration slope and MAPE
# Sample size=308; Prevalence=0.2; C-statistic=0.85; Number of predictors=10
# Calculation takes a few seconds

expected_performance(outcome = "Binary", n = 308, phi = 0.2, c = 0.85, p = 10)

----------------------------  ---------
n                           308.0000
True prevalence               0.2000
True c-statistic              0.8500
Number of predictors         10.0000
---------------------------   0.0000
Mean_CS                       0.8440
SD(CS)                        0.1265
Pr(0.85<CS<1.15)              0.4520
Mean_MAPE                     0.0513
SD(MAPE)                      0.0119
Mean_AUC                      0.8350
SD(AUC)                       0.0092
----------------------------  ---------
```
![image](https://github.com/user-attachments/assets/b334b848-ec07-4fa9-a718-19a355372d11)

As expected, the mean calibration slope for n$rvs=308 is 0.844, smaller than 0.9. The variability is high and translates to 
45% chance of actually getting a model with calibration slope $\in(0.85,1.15)$ when we develop a model with data of that size. Hence, larger size is required.  
In this case, to get a mean calibration slope of 0.9 we need to inflate n$rvs size by 60%! We can confirm that with a sample size of 500 we 
get the desired expected calibration slope:  

``` r
expected_performance(outcome = "Binary", n = 500, phi = 0.2, c = 0.85, p = 10)

----------------------------  ---------
n                           500.0000
True prevalence               0.2000
True c-statistic              0.8500
Number of predictors         10.0000
---------------------------   0.0000
Mean_CS                       0.9020
SD(CS)                        0.0998
Pr(0.85<CS<1.15)              0.6740
Mean_MAPE                     0.0396
SD(MAPE)                      0.0088
Mean_AUC                      0.8410
SD(AUC)                       0.0068
```

![image](https://github.com/user-attachments/assets/d02cda94-a1b5-4618-883d-9e1ed41ec801)

N.B. Although the mean calibration slope is now indeed 0.9 bare in mind that still there is variability in the CS
and *we are not guaranteed* to achieve that performance for every development sample of size 500 ... Indeed, $P(calibration \space slope \in (0.85,1.15))$=67%...


