
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lmtp

<!-- badges: start -->

[![Build
Status](https://travis-ci.com/nt-williams/lmtp.svg?token=DA4a53nWMx6q9LisKdRD&branch=master)](https://travis-ci.com/nt-williams/lmtp)
[![codecov](https://codecov.io/gh/nt-williams/lmtp/branch/master/graph/badge.svg?token=TFQNTischL)](https://codecov.io/gh/nt-williams/lmtp)
[![MIT
license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)
<!-- badges: end -->

> Cross-Validated, Non-parametric Causal Effects Based on Longitudinal
> Modified Treatment Policies

[Nick Williams](https://nicholastwilliams.com) and [Ivan
Diaz](https://idiaz.xyz)

-----

## Scope

`lmtp` is an R package that provides an estimation framework for the
casual effects of longitudinal modified treatment policies described in
\[INSERT IVAN PAPER\]. Two primary estimators are supported, a Targeted
Maximum Likelihood (TML) estimator and a Sequentially Doubly Robust
(SDR) estimator (a substitution and an IPW estimator are provided for
the sake of being thorough but their use is recommended against in favor
of the TML and SDR estimators). In addition to longitudinal effects,
point treatments are naturally supported and both binary and continuous
outcome (both with censoring) are allowed. `lmtp` is built atop the
[`sl3`](https://github.com/tlverse/sl3) package to utilize ensemble
machine learning for estimation.

## Installation

`lmtp` can be installed from GitHub with:

``` r
devtools::install_github("nt-williams/lmtp")
```

For an in-depth look at the package’s functionality, please consult the
accompying
[vignette](https://htmlpreview.github.io/?https://github.com/nt-williams/lmtp/blob/master/vignettes/intro-lmtp.html).

## Example

``` r
library(lmtp)
#> lmtp: Causal Effects Based on Longitudinal Modified Treatment Policies
#> Version: 0.0.7.9000
#> 
library(sl3)
library(future)

# the data: 4 treatment nodes with time varying covariates and a binary outcome
head(sim_t4)
#>   ID L_1 A_1 L_2 A_2 L_3 A_3 L_4 A_4 Y
#> 1  1   2   3   0   1   1   1   1   3 0
#> 2  2   2   1   1   4   0   3   1   2 0
#> 3  3   1   0   1   3   1   2   1   1 1
#> 4  4   1   0   0   3   1   3   1   2 0
#> 5  5   3   3   1   1   0   1   1   2 0
#> 6  6   1   0   0   2   0   3   1   4 0
```

We’re interested in a treatment policy, `d`, where exposure is decreased
by 1 only among subjects whose exposure won’t go below 1 if intervened
upon. The true population outcome under this policy is about 0.305.

``` r
# our treatment policy function to be applied at all time points
d <- function(a) {
  (a - 1) * (a - 1 >= 1) + a * (a - 1 < 1)
}
```

In addition to specifying a treatment policy, we need to specify our
treatment variables, time-varying covariates, and the `sl3` learners to
be used in estimation.

``` r
# our treatment nodes, a character vector of length 4
a <- c("A_1", "A_2", "A_3", "A_4")
# our time varying nodes, a list of length 4
time_varying <- list(c("L_1"), c("L_2"), c("L_3"), c("L_4"))
# our sl3 learner stack: the mean, GLM, and random forest
lrnrs <- make_learner_stack(Lrnr_mean, 
                            Lrnr_glm, 
                            Lrnr_ranger)
```

We’re now ready to estimate the effect of our treatment policy, `d`. In
this example, we’ll use the cross-validated TML estimator with 10 folds.
To speed up computation, we can use parallel processing supported by the
`future` package.

``` r
plan(multiprocess)

lmtp_tmle(sim_t4, a, "Y", time_varying, k = 1, shift = d, 
          learners_outcome = lrnrs, learners_trt = lrnrs, folds = 10)
#> LMTP Estimator: TMLE
#>    Trt. Policy: (d)
#> 
#> Population intervention effect
#>       Estimate: 0.3089
#>     Std. error: 0.0115
#>         95% CI: (0.2863, 0.3315)
```

## Features

| Feature                         |   Status    |
| ------------------------------- | :---------: |
| Point treatment                 |      ✓      |
| Longitudinal treatment          |      ✓      |
| Modified treatment intervention |      ✓      |
| Static intervention             |      ✓      |
| Dynamic intervention            | In progress |
| Continuous treatment            |      ✓      |
| Binary treatment                |      ✓      |
| Continuous outcome              |      ✓      |
| Binary outcome                  |      ✓      |
| Censored outcome                |      ✓      |
| Super learner                   |      ✓      |
| Clustered data                  | In progress |
| Parallel processing             |      ✓      |
| Progress bars                   |      ✓      |

### References
