
#' Simulated longitudinal data
#'
#' A dataset with a binary outcome, four time varying treatment nodes,
#' and four time varying covariates.
#'
#' @format A data frame with 5000 rows and 10 variables:
#' \describe{
#'   \item{ID}{observation ID}
#'   \item{L_1}{Time varying covariate time 1}
#'   \item{A_1}{Treatment node at time 1, effected by L_1}
#'   \item{L_2}{Time varying covariate time 1, effected by L_1 and A_1}
#'   \item{A_2}{Treatment node at time 2, effected by L_2 and A_1}
#'   \item{L_3}{Time varying covariate time 1, effected by L_2 and A_2}
#'   \item{A_3}{Treatment node at time 3, effected by L_3 and A_2}
#'   \item{L_4}{Time varying covariate time 1, effected by L_3 and A_3}
#'   \item{A_4}{Treatment node at time 3, effected by L_4 and A_3}
#'   \item{Y}{Binary outcome at time 5, effected by L_4 and A_4}
#' }
"sim_t4"

#' Simulated longitudinal data with censoring
#'
#' A dataset with a binary outcome, two time varying treatment nodes,
#' two time varying covariates, and two censoring indicators.
#'
#' @format A data frame with 1000 rows and 10 variables:
#' \describe{
#'   \item{L1}{Time varying covariate time 1}
#'   \item{A1}{Treatment node at time 1, effected by L_1}
#'   \item{C1}{Censoring indicator that the observation is observed after time 1}
#'   \item{L2}{Time varying covariate at time 2, effected by L_1 and A_1}
#'   \item{A2}{Treatment node at time 2, effected by L_2 and A_1}
#'   \item{C2}{Censoring indicator that the observation is observed after time 2}
#'   \item{Y}{Binary outcome at time 3, effected by L_2 and A_2}
#' }
"sim_cens"
