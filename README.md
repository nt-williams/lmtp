
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lmtp

<!-- badges: start -->

[![Build
Status](https://travis-ci.com/nt-williams/lmtp.svg?token=DA4a53nWMx6q9LisKdRD&branch=master)](https://travis-ci.com/nt-williams/lmtp)
[![codecov](https://codecov.io/gh/nt-williams/lmtp/branch/master/graph/badge.svg?token=TFQNTischL)](https://codecov.io/gh/nt-williams/lmtp)
<!-- badges: end -->

> Non-parametric Causal Effects Based on Longitudinal Modified Treatment
> Policies

# Installation

The development version can be installed from GitHub with:

``` r
devtools::install_github("nt-williams/lmtp")
```

# Example

Simulate some data…

``` r
set.seed(429)
n_obs <- 1000
L1 <- rbinom(n_obs, 1, 0.5)
A1 <- rnorm(n_obs, mean = 2 * L1, sd = 1)
L2 <- rbinom(n_obs, 1, plogis(A1 + L1 + rnorm(n_obs, mean = 0, sd = 1)))
A2 <- rnorm(n_obs, mean = 2 * A1 + L2, sd = 1)
Y <- rbinom(n_obs, 1, plogis(A2 + L2 + rnorm(n_obs, mean = 0, sd = 1)))
df <- data.frame(L1, A1, L2, A2, Y)
a <- c("A1", "A2")
nodes <- list(c("L1"),
              c("L2"))

head(df)
#>   L1        A1 L2        A2 Y
#> 1  1 2.2500683  1 4.9955517 1
#> 2  0 1.7413720  1 5.1433794 1
#> 3  1 0.5989326  0 0.8569932 1
#> 4  1 3.5842336  1 7.5588927 1
#> 5  0 0.9392441  1 2.5614636 1
#> 6  0 1.4205120  1 2.4504096 1
```

An intervention where exposure is increased by 0.5 at all time points
for all subjects. The true value under this intervetion is about 0.88

``` r
library(lmtp)
#> lmtp: Causal Effects Based on Longitudinal Modified Treatment Policies
#> Version: 0.0.1.9000

lmtp_tmle(df, a, "Y", nodes, k = 0, shift = function(x) x + 0.5,
          outcome_type = "binomial",
          learner_stack_Q = sl3::make_learner(sl3::Lrnr_glm),
          learner_stack_g = sl3::make_learner(sl3::Lrnr_glm))
#>   Estimating propensity [=====================================] 100%  0s
#>   Estimating regression [=====================================] 100%  0s
#> 
#> Estimator: TMLE
#> 
#> Population intervention effect
#>    Estimate: 0.8854
#>  Std. error: 0.0098
#>      95% CI: (0.8662, 0.9046)
```
