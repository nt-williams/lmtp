#' Simulated Longitudinal Data
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

#' Simulated Longitudinal Data With Censoring
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

#' Simulated Point-treatment Survival Data
#'
#' A dataset with a time-to-event outcome, two baseline nodes, a binary
#' point treatment, six past-time outcome nodes, and six censoring indicators.
#'
#' @format A data frame with 2000 rows and 16 variables:
#' \describe{
#'   \item{W1}{Binary baseline variable.}
#'   \item{W2}{Categorical baseline variable.}
#'   \item{trt}{Binary treatment variable.}
#'   \item{C.0}{Censoring indicator that the observation is observed future time points.}
#'   \item{Y.1}{Outcome node at time 1.}
#'   \item{C.1}{Censoring indicator that the observation is observed future time points.}
#'   \item{Y.2}{Outcome node at time 2.}
#'   \item{C.2}{Censoring indicator that the observation is observed future time points.}
#'   \item{Y.3}{Outcome node at time 3.}
#'   \item{C.3}{Censoring indicator that the observation is observed future time points.}
#'   \item{Y.4}{Outcome node at time 4.}
#'   \item{C.4}{Censoring indicator that the observation is observed future time points.}
#'   \item{Y.5}{Outcome node at time 5.}
#'   \item{C.5}{Censoring indicator that the observation is observed future time points.}
#'   \item{Y.6}{Final outcome node.}
#' }
"sim_point_surv"

#' Simulated Time-varying Survival Data
#'
#' A dataset with a time-to-event outcome, one baseline nodes, two time-varying
#' covariates, a binary time-varying treatment, two outcome nodes,
#' and two censoring indicators. Data-generating mechanism taken from
#' Lendle, Schwab, Petersen, and van der Laan (<https://www.jstatsoft.org/article/view/v081i01>).
#'
#' @format A data frame with 500 rows and 11 variables:
#' \describe{
#'   \item{L0.a}{Continuous baseline variable.}
#'   \item{L0.b}{Time varying covariate at baseline.}
#'   \item{L0.c}{Time varying covariate at baseline.}
#'   \item{A0}{Treatment variable at baseline}
#'   \item{C0}{Censoring indicator that the observation is observed future time points.}
#'   \item{L1.a}{Time varying covariate at time 1.}
#'   \item{L1.b}{Time varying covariate at time 1.}
#'   \item{Y1}{Outcome node at time 1.}
#'   \item{A1}{Treatment variable at time 1.}
#'   \item{C1}{Censoring indicator that the observation is observed future time points.}
#'   \item{Y2}{Final outcome node.}
#' }
"sim_timevary_surv"

#' Simulated Multivariate Exposure Data
#'
#' A dataset with a continuous outcome, three baseline covariates,
#' and two treatment variables.
#'
#' @format A data frame with 2000 rows and 6 variables:
#' \describe{
#'   \item{C1}{Continuous baseline variable.}
#'   \item{C2}{Continuous baseline variable.}
#'   \item{C3}{Continuous baseline variable.}
#'   \item{D1}{Treatment variable one at baseline.}
#'   \item{D2}{Treatment variable two at baseline.}
#'   \item{Y}{Continuous outcome}
#' }
"multivariate_data"
