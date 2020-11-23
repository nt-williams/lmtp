# lmtp 1.0.0.9000

* Now relies on SuperLearner for estimation
* fixed a bug that outputted outcome regressions and density ratios in incorrect order compared to the original data
* fixed a bug in the missing data check that threw an error for missing data 
 after an observation experiences the outcome
* fixed a bug that resulted in NA censoring indicators throwing an error for missing data
* Added a `NEWS.md` file to track changes to the package.
* license change to GPL-3
