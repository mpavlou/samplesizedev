
<!-- README.md is generated from README.Rmd. Please edit that file -->

# samplesizedev : Simulation-based sample size calculations for the development of risk models 

<!-- badges: start -->
<!-- badges: end -->


R package to calculate the sample size for the development of risk models for binary outcomes.

Last Update: 17/06/2026

Related papers: 

**"An evaluation of sample size requirements for developing risk prediction models with binary outcomes"**, BMC Medical Research Methodology https://doi.org/10.1186/s12874-024-02268-5

$\textcolor{#f00}{\large  \textbf{NEW preprint:}}$ **"Sample Size Calculations for the Development of Risk Prediction Models that Account for Performance Variability"**,
https://doi.org/10.48550/arXiv.2509.14028

**"Penalized Regression Methods With Modified Cross-Validation and Bootstrap Tuning Produce Better Prediction Models"**. Biometrical Journal, https://onlinelibrary.wiley.com/doi/10.1002/bimj.202300245

### Why do we need 'samplesizedev'?

Riley et al. (2019) proposed calculating the development for a prediction model to ensure that the expected overfitting, as quantified by the calibration slope, is small . The default target is calibration slope=0.9. Other criteria which include the optimism in R2 Nagelkerke and the precision in the mean predicted risk. The formula which aims to control model overfitting ('calibration' formula - C1) most often gives that highest sample size and our article we focused primarily around this formula. While the calibration formula performed well for models with C-statistic/C-index<0.8, we found that that it substantially ***underestimated*** the sample size when the predictive strength of the model was higher. The sample sizes often needed to be increased by 50% or even doubled to meet the calibration targets.

Hence, we developed the **new package 'samplesizedev'** which performs ***unbiased sample size calculations*** regardless of model strength. Our software uses simulation in the background so calculations can take a bit to run (from  few seconds to 3 minutes depending on the specific scenario and computational power). Currently it can be used for the development of risk models for binary outcomes; functionality for ***time to event outcomes***  is under development. 

### $\textcolor{#f00}{\large  \textbf{UPDATES}}$

The package has been recently been updated (April 2026) to:
- perform sample size calculations than also ***control the*** $\textcolor{#f00}{variability}$ ***in the calibration slope, instead of *just* the expected value***.
- perform sample size calculations using analytical expressions which result in a very good approximation and impressive speed improvements 
- obtain the sampling distribution for a variety of performance measures (e.g. C-statistic, Brier score etc) for a given sample size
- obtain the sampling distribution for individual predicted probabilities for a given sample size  

**The importance of accounting for variability in performance.** The first above is very important because, as shown in our recent [preprint on performance variability](https://doi.org/10.48550/arXiv.2509.14028), under the current sample size calculations aiming at expected calibration slope (CS) of 0.9,  the variability in performance can be very high when the number of predictors is small. Therefore, while it might be thought that a small model developed with a relatively small sample size can avoid model overfitting, this can be misleading. That's because even though performance is controlled on average,  e.g., expected calibration slope E(CS) = 0.9, the variability can be very high. In this context, E(CS)=0.9 is interpreted to mean that if one were to collect many datasets of the recommended size, and  validate them on large external dataset, then the calibration slope would be on average around 0.9. However, if the variability is very high, the probability of actually obtaining an individual dataset with calibration slope close to 0.9 might be unacceptably low (see examples in the paper and below). 

Hence, in our more recent work we $\textcolor{#f00}{\text{develop sample size calculations where we aim to control the probability of acceptable performance, rather}}$ $\textcolor{#f00}{ \text{than just performance on average}}$. In a simulation-based framework this approach can be easily implemented for any performance metric and a suitably defined range of acceptable performance. Here and in the paper above, we focused on the calibration slope. In addition to the simulation-based approach, we also derived an approximate analytical calculation that avoids simulation and is, hence,  very quick. 


### How does 'samplesizedev' work?

The software requires information on the anticipated values of the:
- outcome prevalence
- c-statistic (AUC)
- number of predictor variables

Based on the characteristics above it can perform actions based on two core functions:

1. **Calculate the required sample size** to achieve a target expected calibration slope or Mean Absolute Prediction Error (MAPE) (function **'samplesizedev'**)

2. $\textcolor{#f00}{\textbf{NEW:}}$ **Calculate the required sample size** to achieve a $\textcolor{#f00}{\text{a high  probability  of  a   model   with   acceptable   calibration}}$ (function **'samplesizedev'**)

3. **Calculate the expected model performance** at a given sample size (function **'expected_performance'**)

4. $\textcolor{#f00}{\text{NEW:}}$ the  **sample size** for E(CS) or PrAP(CS) using a quick implementation (no simulation) based on analytical approximations (option 'quick = TRUE')



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

1.  the **sample size** to achieve a target expected calibration slope (e.g. target expected calibration slope **E(CS)=0.9**)  or 

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

samplesizedev(outcome = "Binary", S = 0.9, phi = 0.2, c = 0.85, p = 10, quick = FALSE)
> [1] "Optimisation Started: check progress on the appearing plots..."
> $rvs
> [1] 309
> $analytical_corrected
> [1] 486
> $sim
> [1] 535
$note
[1] "Monte Carlo Simulation Error (MCSE) = 0.003. The sample size within 1 MCSE from S=0.9 would be approximately 517-553. For lower MCE error increase the number of simulations."

# $sim is the sample size calculated by simulation
# $rvs is the sample size calculated using the approach of Riley et al. (2019) (RvS formula Criterion 1 - overfitting)
# $analytical_corrected is the sample size after applying bias reduction to Riley's calibration formula
# Note: Check the Monte Carlo simulation error on the appearing plots

```

Alternatively, one may use the 'quick = TRUE' option which directly uses the bias-reduction adjustment for Riley's formula, which uses the 'adjusted C-statistic' as input:

``` r
# Fast calculation (approximation) 

samplesizedev(outcome = "Binary", S = 0.9, phi = 0.2, c = 0.85, p = 10, quick = TRUE)
> $rvs
> [1] 309
> $analytical_corrected
> [1] 486
``` 
We note that the adjusted C-statistic value used for the corrected calculation above is 0.79:

``` r
c_adj(target.prev = 0.2, target.c = 0.85, p = 10)
> [1] 0.788
```

The sample size calculated using simulation is n$sim=535 which corresponds to CS=0.9. In comparison, 
the sample size using previously proposed formulae is n$rvs=309. According to the findings in our paper
the RvS overfitting formula underestimates the sample size for high C-statistic. Thus, the expected calibration slope will
be in fact lower than we had aim for this size.  We can verify this using the second command of our package, 'expected_performance'.


#### Calculation of expected model performance (CS, C-statistic, MAPE etc and their variability) for a given sample size and model characteristics, and variability in individual predicted probabilities

``` r
# Calculate the expected calibration slope and MAPE
# Sample size=308; Prevalence=0.2; C-statistic=0.85; Number of predictors=10
# Calculation takes a few seconds

expected_performance(outcome = "Binary", n = 309, phi = 0.2, c = 0.85, p = 10, method = "MLE")

n                                     309.0000
True prevalence                         0.2000
True c-statistic                        0.8500
Number of predictors                   10.0000
-------------------------------------  
Sampling distribution of               
aggregate performance measures          
-------------------------------------   
Mean Calibration Slope (CS)             0.8370
SD(CS)                                  0.1229
Pr(0.85<CS<1.15)                        0.4300
Mean MAPE                               0.0509
SD(MAPE)                                0.0117
Mean_AUC                                0.8360
SD(AUC)                                 0.0085
SD(Average Predicted Risk)              0.0190
Median CS                               0.8310
Mean Brier Score                        0.1180
Sensitivity (threshold prevalence)      0.7860
Net Benefit (threshold prevalence)      0.1120
------------------------------------    
Sampling distribution of Individual     
Predicted Probability (IPP) = 0.112     
------------------------------------    
Median IPP                              0.1020
SD(IPP)                                 0.0630
------------------------------------
```
<img width="788" height="586" alt="github_fig1" src="https://github.com/user-attachments/assets/82ed1f9a-eb00-4c0a-a81d-4a3c2308f250" />


As expected, the median calibration slope for n$rvs=309 is 0.831, smaller than 0.9. The variability is high and translates to 
43% chance of actually getting a model with calibration slope $\in(0.85,1.15)$ when we develop a model with data of that size. Hence, larger size is required.  
In this case, to get a median calibration slope of 0.9 we need to inflate n$rvs size by approximately 60%! We can confirm that with a sample size of 535 we 
get the desired median calibration slope:  

``` r
expected_performance(outcome = "Binary", n = 535, phi = 0.2, c = 0.85, p = 10, method = "MLE")

n                                     535.0000
True prevalence                         0.2000
True c-statistic                        0.8500
Number of predictors                   10.0000
-------------------------------------  
Sampling distribution of                
aggregate performance measures          
-------------------------------------   
Mean Calibration Slope (CS)             0.9070
SD(CS)                                  0.0962
Pr(0.85<CS<1.15)                        0.7100
Mean MAPE                               0.0383
SD(MAPE)                                0.0081
Mean AUC                                0.8420
SD(AUC)                                 0.0057
SD(Average Predicted Risk)              0.0150
Median CS                               0.9000
Median Brier Score                      0.1150
Median Sensitivity (threshold prev)     0.8040
Median Net Benefit (threshold prev)     0.1150
------------------------------------    
Sampling distribution of Individual     
Predicted Probability (IPP) = 0.112     
------------------------------------   
Median IPP                              0.1090
SD(IPP)                                 0.0440
------------------------------------    
```

<img width="788" height="586" alt="github_fig2" src="https://github.com/user-attachments/assets/6127960f-e82b-4719-a623-d80cf172fcd8" />



N.B. Although the mean calibration slope is now indeed 0.9 bare in mind that still there is variability in the CS
and *we are not guaranteed* to achieve that performance for every development sample of size 500 ... Indeed, $P(calibration \space slope \in (0.85,1.15))$=67%...

 $\textcolor{#f00}{\textbf{NEW:}}$ **Alternative fitting methods.** The package now provides the option to calculation expected performance after using shrinkage with a) post estimation shrinkage via the linear shrinkage factor (method="LSF"), b) modified Ridge (method = "ridge") and c) modified LASSO (method = "lasso"), where  [modified Ridge and LASSO were proposed by Pavlou et al. (2024)](https://onlinelibrary.wiley.com/doi/10.1002/bimj.202300245)  to solve problems with increase variability of standard Ridge and LASSO. Note that the implementation for these alternative methods is time consuming as the corresponding algorithms are much slower than MLE (about 5-10 times slower).
 
``` r
# Try the alternative estimation methods ridge and LASSO
expected_performance(outcome = "Binary", n = 535, phi = 0.2, c = 0.85, p = 10, method = "ridge")
                                          [,1]
n                                     535.0000
True prevalence                         0.2000
True c-statistic                        0.8500
Number of predictors                   10.0000
-------------------------------------  
Sampling distribution of                
aggregate performance measures          
-------------------------------------   
Mean Calibration Slope (CS)             1.0160
SD(CS)                                  0.1123
Pr(0.85<CS<1.15)                        0.8200
Mean MAPE                               0.0375
SD(MAPE)                                0.0082
Mean AUC                                0.8420
SD(AUC)                                 0.0055
SD(Average Predicted Risk)              0.0150
Median CS                               1.0070
Median Brier Score                      0.1150
Median Sensitivity (threshold prev)     0.7690
Median Net Benefit (threshold prev)     0.1040
------------------------------------    
Sampling distribution of Individual     
Predicted Probability (IPP) = 0.112     
------------------------------------    
Median IPP                              0.1200
SD(IPP)                                 0.0430
------------------------------------    

```


#### Calculation of sample size for given model characteristics, aiming at $\textcolor{#f00}{ \textbf{Probability of acceptable calibration, PrAP(CS)=0.8}}$ 

``` r
# Calculate the sample size Size for Probability of Acceptable Performance (PrAP=0.8),
# where Acceptable Performance is defined $S\in (0.85, 1.15)$
# Performance target: PrAP(S)=0.8; Prevalence=0.2; c-statistic=0.85; Number of predictors=10

samplesizedev(outcome="Binary", l_s= 0.85, u_s = 1.15, PrAP_s = 0.8,
              phi = 0.2, c = 0.85, p = 10,  quick = FALSE)

$sim
[1] 701
$analytical
[1] 649

# $sim is the sample size calculated by simulation to ensure that PrAP(S)=0.8
# $analytical is the sample size calculated using an analytical approximation to ensure that PrAP(S)=0.8

```
The sample size calculated using simulation targeting at E(S)=0.9 is 535, while the sample size to ensure that PrAP(S)=0.8 is 701.

As before, one may use the 'quick = TRUE' option which uses a bias-reduction analytical method. This runs much faster and provides a very good approximation:
``` r
samplesizedev(outcome="Binary", l_s= 0.85, u_s = 1.15, PrAP_s = 0.8,
              phi = 0.2, c = 0.85, p = 10,  quick = TRUE)

$analytical
[1] 649
```



