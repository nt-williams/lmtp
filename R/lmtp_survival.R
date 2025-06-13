#' LMTP Survival Curve Estimator
#'
#' Wrapper around \code{lmtp_tmle} and \code{lmtp_sdr} for survival outcomes to estimate the entire survival curve.
#' Estimates are reconstructed using isotonic regression to enforce monotonicity of the survival curve.
#' \bold{Confidence intervals correspond to marginal confidence intervals for the survival curve, not simultaneous intervals.}
#'
#' @param data \[\code{data.frame}\]\cr
#'  A \code{data.frame} in wide format containing all necessary variables
#'  for the estimation problem. Must not be a \code{data.table}.
#' @param trt \[\code{character}\] or \[\code{list}\]\cr
#'  A vector containing the column names of treatment variables ordered by time.
#'  Or, a list of vectors, the same length as the number of time points of observation.
#'  Vectors should contain column names for the treatment variables at each time point. The list
#'  should be ordered following the time ordering of the model.
#' @param outcomes \[\code{character}\]\cr
#'  A vector containing the columns names of intermediate outcome variables and the final
#'  outcome variable ordered by time. Only numeric values are allowed. Variables should be coded as 0 and 1.
#' @param baseline \[\code{character}\]\cr
#'  An optional vector containing the column names of baseline covariates to be
#'  included for adjustment at every time point.
#' @param time_vary \[\code{list}\]\cr
#'  A list the same length as the number of time points of observation with
#'  the column names for new time-varying covariates introduced at each time point. The list
#'  should be ordered following the time ordering of the model.
#' @param cens \[\code{character}\]\cr
#'  An optional vector of column names of censoring indicators the same
#'  length as the number of time points of observation. If missingness in the outcome is
#'  present or if time-to-event outcome, must be provided.
#' @param compete \[\code{character}\]\cr
#'  An optional vector of column names of competing risk indicators the same
#'  length as the number of time points of observation.
#' @param shift \[\code{closure}\]\cr
#'  A two argument function that specifies how treatment variables should be shifted.
#'  See examples for how to specify shift functions for continuous, binary, and categorical exposures.
#' @param shifted \[\code{data.frame}\]\cr
#'  An optional data frame, the same as in \code{data}, but modified according
#'  to the treatment policy of interest. If specified, \code{shift} is ignored.
#' @param estimator \[\code{character(1)}\]\cr
#'  The estimator to use. Either \code{"lmtp_tmle"} or \code{"lmtp_sdr"}.
#' @param k \[\code{integer(1)}\]\cr
#'  An integer specifying how previous time points should be
#'  used for estimation at the given time point. Default is \code{Inf},
#'  all time points.
#' @param mtp \[\code{logical(1)}\]\cr
#'  Is the intervention of interest a modified treatment policy?
#'  Default is \code{TRUE}. If treatment variables are continuous this should be \code{TRUE}.
#'  If the treatment variables are not continuous, setting this to \code{FALSE} may improve performance.
#' @param id \[\code{character(1)}\]\cr
#'  An optional column name containing cluster level identifiers.
#' @param learners_outcome \[\code{character}\]\cr A vector of \code{SuperLearner} algorithms for estimation
#'  of the outcome regression. Default is \code{"SL.glm"}, a main effects GLM.
#' @param learners_trt \[\code{character}\]\cr A vector of \code{SuperLearner} algorithms for estimation
#'  of the exposure mechanism. Default is \code{"SL.glm"}, a main effects GLM.
#'  \bold{Only include candidate learners capable of binary classification}.
#' @param folds \[\code{integer(1)}\]\cr
#'  The number of folds to be used for cross-fitting.
#' @param weights \[\code{numeric(nrow(data))}\]\cr
#'  An optional vector containing sampling weights.
#' @param control \[\code{list()}\]\cr
#'  Output of \code{lmtp_control()}.
#'
#' @return A list of class \code{lmtp_survival} containing \code{lmtp} objects for each time point.
#'
#' @example inst/examples/lmtp_survival-ex.R
#' @export
lmtp_survival <- function(data, trt, outcomes, baseline = NULL, time_vary = NULL,
                          cens = NULL, compete = NULL,
                          shift = NULL, shifted = NULL,
                          estimator = c("lmtp_tmle", "lmtp_sdr"),
                          k = Inf, mtp = TRUE, id = NULL,
                          learners_outcome = "SL.glm",
                          learners_trt = "SL.glm",
                          folds = 10,
                          weights = NULL,
                          control = lmtp_control()) {

  checkmate::assertCharacter(outcomes, min.len = 2, null.ok = FALSE, unique = TRUE, any.missing = FALSE)

  estimator <- match.arg(estimator)
  time_horizon <- length(outcomes)
  estimates <- vector("list", time_horizon)

  args <- list(
    data = data,
    baseline = baseline,
    shift = shift,
    shifted = shifted,
    k = k,
    mtp = mtp,
    id = id,
    learners_outcome = learners_outcome,
    learners_trt = learners_trt,
    folds = folds,
    weights = weights,
    control = control
  )

  if (length(trt) == 1) args$trt <- trt
  if (length(time_vary) == 1) args$time_vary <- time_vary

  time <- 1
  cli::cli_progress_step("Working on time {time}/{time_horizon}...")
  for (time in seq_len(time_horizon)) {
    if (length(trt) > 1) args$trt <- trt[seq_len(time)]
    if (length(args$time_vary) > 1) args$time_vary <- time_vary[seq_len(time)]
    args$outcome <- outcomes[seq_len(time)]
    args$cens <- cens[seq_len(time)]
    args$compete <- compete[seq_len(time)]
    args$outcome_type <- ifelse(time == 1, "binomial", "survival")

    if (estimator == "lmtp_tmle") {
      estimates[[time]] <- future::future(do.call(lmtp_tmle, args), seed = TRUE)
    } else {
      estimates[[time]] <- future::future(do.call(lmtp_sdr, args), seed = TRUE)
    }

    cli::cli_progress_update()
  }

  cli::cli_progress_done()
  estimates <- future::value(estimates)
  estimates <- fix_surv_time1(estimates)
  estimates <- isotonic_projection(estimates)

  class(estimates) <- "lmtp_survival"
  estimates
}

isotonic_projection <- function(x) {
  estim <- do.call("rbind", lapply(x, tidy))
  iso_fit <- isotone::gpava(1:length(x), 1 - estim$estimate)
  for (i in seq_along(x)) {
    x[[i]]$estimate <- ife::ife(1 - iso_fit$x[i],
                                x[[i]]$estimate@eif,
                                x[[i]]$estimate@weights,
                                x[[i]]$estimate@id)
  }
  x
}
