
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

We have a data set with treatment nodes at 4 time points and a binary
outcome at time 5.

``` r
head(lmtp::sim_t4)
#>   ID L_1 A_1 L_2 A_2 L_3 A_3 L_4 A_4 Y
#> 1  1   2   1   0   1   1   3   1   2 0
#> 2  2   2   3   1   3   1   3   1   4 1
#> 3  3   3   2   1   2   1   2   0   3 0
#> 4  4   1   0   0   4   1   1   1   4 1
#> 5  5   3   2   1   3   1   3   1   3 0
#> 6  6   3   4   0   2   1   3   1   0 0
```

We’re interested in the effect of an intervention on outcome `Y` where
treatment is decreased by 1 at all time points for observations whose
treatment value won’t go below 0 if intervened upon. The true value
under this intervention is about 0.48.

![](https://gist.githubusercontent.com/nt-williams/2488fef9e94c7ef1a3920c2682433980/raw/2a339c10d651aa18dc48188b4017c605c30e2405/lmtp-readme-example.svg?sanitize=true)
