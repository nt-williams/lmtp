# lmtp 1.0.0.9000

* Now relies on the **SuperLearner** package for estimation.
* Survival outcome now estimates survival probability instead of cumulative incidence. This fixes a bug with IPW and survival problems.
* Outcome type now accepts "survival" for explicit indication of a survival outcome
* New `.trimming` parameter for trimming extreme density ratios.
* Fixed a bug that outputted outcome regressions and density ratios in incorrect order compared to the original data.
* Fixed a bug in the missing data check that threw an error for missing data after an observation experiences the outcome.
* Fixed a bug that resulted in `NA` censoring indicators throwing an error for missing data.
* Added a `NEWS.md` file to track changes to the package.
* license change to GPL-3
