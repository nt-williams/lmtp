
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lmtp <img src='man/figures/lmtp.png' align="right" height="139" /></a>

<!-- badges: start -->

<!-- [![Build Status](https://travis-ci.com/nt-williams/lmtp.svg?token=DA4a53nWMx6q9LisKdRD&branch=master)](https://travis-ci.com/nt-williams/lmtp) -->

<!-- [![codecov](https://codecov.io/gh/nt-williams/lmtp/branch/master/graph/badge.svg)](https://codecov.io/gh/nt-williams/lmtp) -->

<!-- [![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT) -->

<!-- [![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active) -->

<!-- [![DOI](https://zenodo.org/badge/251356023.svg)](https://zenodo.org/badge/latestdoi/251356023) -->

<!-- badges: end -->

> Non-parametric Causal Effects of Feasible Interventions Based on
> Modified Treatment Policies

Nick Williams and Ivan Diaz

-----

## Installation

`lmtp` can be installed from GitHub with:

``` r
devtools::install_github("nt-williams/lmtp")
```

The stable, development version can be installed from GitHub with:

``` r
devtools::install_github("nt-williams/lmtp@devel")
```

## Scope

`lmtp` is an R package that provides an estimation framework for the
casual effects of feasible interventions based on point-treatment and
longitudinal modified treatment policies as described in Diaz, Williams,
Hoffman, and Schenck (2020). Two primary estimators are supported, a
targeted maximum likelihood (TML) estimator and a sequentially doubly
robust (SDR) estimator (a G-computation and an inverse probability of
treatment weighting estimator are provided for the sake of being
thorough but their use is recommended against in favor of the TML and
SDR estimators). Both binary and continuous outcomes (both with
censoring) are allowed. `lmtp` is built atop the
[`sl3`](https://github.com/tlverse/sl3) package to utilize ensemble
machine learning for estimation. The treatment mechanism is estimated
via a density ratio classification procedure irrespective of treatment
variable type providing decreased computation time when treatment is
continuous. Dynamic treatment regimes are also supported.

For an in-depth look at the package’s functionality, please consult the
accompanying
[article](https://htmlpreview.github.io/?https://gist.githubusercontent.com/nt-williams/ddd44c48390b8d976fad71750e48d8bf/raw/45db700a02bf92e2a55790e60ed48266a97ca4e7/intro-lmtp.html).

### Features

| Feature                         | Status |
| ------------------------------- | :----: |
| Point treatment                 |   ✓    |
| Longitudinal treatment          |   ✓    |
| Modified treatment intervention |   ✓    |
| Static intervention             |   ✓    |
| Dynamic intervention            |   ✓    |
| Continuous treatment            |   ✓    |
| Binary treatment                |   ✓    |
| Categorical treatment           |   ✓    |
| Missingness in treatment        |        |
| Continuous outcome              |   ✓    |
| Binary outcome                  |   ✓    |
| Censored outcome                |   ✓    |
| Mediation                       |        |
| Super learner                   |   ✓    |
| Clustered data                  |   ✓    |
| Parallel processing             |   ✓    |
| Progress bars                   |   ✓    |

## Example

``` r
library(lmtp)

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
policy <- function(data, trt) {
  (data[[trt]] - 1) * (data[[trt]] - 1 >= 1) + data[[trt]] * (data[[trt]] - 1 < 1)
}
```

In addition to specifying a treatment policy, we need to specify our
treatment variables and time-varying covariates.

``` r
# our treatment nodes, a character vector of length 4
a <- c("A_1", "A_2", "A_3", "A_4")
# our time varying nodes, a list of length 4
time_varying <- list(c("L_1"), c("L_2"), c("L_3"), c("L_4"))
```

We can now estimate the effect of our treatment policy, `d`. In this
example, we’ll use the cross-validated TML estimator with 10 folds.

``` r
lmtp_tmle(sim_t4, a, "Y", time_vary = time_varying, k = 0, shift = policy, folds = 10)
#> LMTP Estimator: TMLE
#>    Trt. Policy: (policy)
#> 
#> Population intervention effect
#>       Estimate: 0.2598
#>     Std. error: 0.0188
#>         95% CI: (0.2228, 0.2967)
```

## Similar Implementations

A variety of other R packages perform similar tasks as `lmtp`. However,
`lmtp` is the only R package currently capable of estimating causal
effects for binary, categorical, and continuous exposures in both the
point treatment and longitudinal setting using traditional causal
effects or modified treatment policies.

  - [`txshift`](https://github.com/nhejazi/txshift)  
  - [`tmle3`](https://github.com/tlverse/tmle3)  
  - [`tmle3shift`](https://github.com/tlverse/tmle3shift)
  - [`ltmle`](https://CRAN.R-project.org/package=ltmle)  
  - [`tmle`](https://CRAN.R-project.org/package=tmle)

## Citation

Please cite the following when using `lmtp` in publications. Citation
should include both the R package and the paper establishing the
statistical methodology.

    @Manual{,
        title = {lmtp: {Non}-parametric {Causal} {Effects} of {Feasible} {Interventions} {Based} on {Modified} {Treatment} {Policies}},
        author = {Nicholas T Williams and Iván Díaz},
        year = {2020},
        note = {R package version 0.0.5},
        doi = {10.5281/zenodo.3874931}, 
        url = {https://github.com/nt-williams/lmtp}
    }
    
    @Article{,
        journal = {arxiv},
        title = {Non-parametric causal effects based on longitudinal modified treatment policies},
        author = {Iván Díaz and Nicholas Williams and Katherine L Hoffman and Edward J Schneck},
        year = {2020},
        url = {https://arxiv.org/abs/2006.01366v2}
    }

## References

Diaz I, Williams N, Hoffman KL, Schenck, EJ (2020). *Non-Parametric
Causal Effects Based on Longitudinal Modified Treatment Policies*.
<https://arxiv.org/abs/2006.01366>
