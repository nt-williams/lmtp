#' @title Flip SDR Estimator
#' @description
#' Estimate the “flip” effect via a sequentially doubly robust (SDR) procedure that 
#' applies either smooth trimming or overlap weighting at each time point.
#'
#' @param data               A \code{data.frame} containing all relevant variables.
#' @param trt                Character vector of treatment column names (time‐varying, each ending in \code{"_tXX"}).
#' @param outcome            Name of the outcome column (character).
#' @param baseline           Character vector of baseline covariate names, or \code{NULL} if none.
#' @param time_vary          Character vector of time‐varying covariate names (each ending in \code{"_tXX"}), or \code{NULL}.
#' @param cens               Name of the censoring indicator column (character), or \code{NULL} if no censoring.
#' @param compete            Name of the competing event indicator column (character), or \code{NULL} if none.
#' @param target_regime      Vector or list specifying the target treatment regime at each time point.
#'                           Must have length equal to \code{length(trt)} (or \code{task$tau} internally).
#' @param overlap            Logical; if \code{TRUE}, use overlap weights, otherwise use smooth trimming.
#' @param trimming_threshold Numeric threshold for smooth trimming. Values outside
#'                           this threshold receive zero weight under smooth trimming.
#'                           Ignored if \code{overlap = TRUE}.
#' @param smoothing_constant Numeric smoothing parameter (larger values --> weight function is closer to a 
#'                           hard indicator). Ignored if \code{overlap = TRUE}.
#' @param k                  Integer (or \code{Inf}); history window size. If \code{Inf}, use the full history.
#' @param outcome_type       Character; one of \code{c("binomial", "continuous", "survival")}.
#'                           Determines how pseudo‐outcomes are constructed and whether bounding applies.
#' @param id                 Name of the subject‐identifier column (character), or \code{NULL}.  
#'                           If \code{NULL}, no ID is stored in the returned \code{ife} object.
#' @param bounds             Optional numeric vector of length 2 (\code{c(lower, upper)}) giving 
#'                           outcome bounds. If \code{NULL}, no bounding is enforced on regression fits.
#' @param learners_outcome   Character vector of SuperLearner (or other) learners to use for the
#'                           sequential outcome regression steps.
#' @param learners_trt       Character vector of SuperLearner (or other) learners to use for
#'                           the time‐varying treatment (propensity) models.
#' @param folds              Integer ≥ 2; number of cross‐fitting folds for both propensity and
#'                           outcome regressions.
#' @param weights            Optional numeric vector of sample weights (length = \code{nrow(data)}).
#'                           If \code{NULL}, no sample weighting is applied.
#' @param control            An \code{lmtp_control} object controlling learner options, 
#'                           parallelization, and other fitting parameters.
#' @param life               Logical; if \code{TRUE}, the returned list will also include
#'                           \code{$task} (the \code{LmtpTask} object) and \code{$control}, so that
#'                           one can continue estimation via \code{flip_sdr_second_PO} or
#'                           \code{flip_sdr_treatment}. If \code{FALSE} (default), only final
#'                           outputs are returned.
#'
#' @return A list of class \code{flip} containing:
#' \describe{
#'   \item{\code{estimator}}{Character; the string \code{"SDR"} indicating this is the SDR flip estimator.}
#'   \item{\code{ifvalues}}{Numeric vector of length \code{nrow(data)} containing the final
#'                          pseudo‐outcome (influence‐function) values at time $t = 1$ after back‐substitution.}
#'   \item{\code{estimate}}{An \code{\link[ife]{ife}} object with 
#'                          \code{x = mean(ifvalues)} and \code{eif = ifvalues}, indexed by \code{id} if provided.}
#'   \item{\code{prop_scores}}{Numeric matrix of dimension $n \times \tau$, containing the cross‐fitted
#'                          estimated propensity scores at each time point.}
#'   \item{\code{fits_prop}}{A list of length \code{folds}, where each entry is itself a list of fitted
#'                          propensity‐model objects (one per time point) from \code{cf_propensity}.}
#'   \item{\code{fits_m}}{A list of length $\tau$, where each entry is the list of fitted outcome‐regression
#'                        objects (one per fold) returned by \code{cf_sequential_regression()} at that time point.}
#'   \item{\code{outcome_type}}{Character; same as the \code{outcome_type} argument, indicating \code{"binomial"}, 
#'                             \code{"continuous"}, or \code{"survival"}.}
#'   \item{\code{task}}{(Only if \code{life = TRUE}.) The \code{LmtpTask} object constructed internally,
#'                      containing original data, shift plan, and relevant metadata (used for continuation).}
#'   \item{\code{control}}{(Only if \code{life = TRUE}.) The \code{lmtp_control} object used in this run.}
#' }
#'
#' @details
#' The \code{flip_sdr} function implements a sequentially doubly robust estimator for a “flip” intervention:
#' for each time $t = 1, \dots, \tau$, one either trims units whose estimated density ratio falls below 
#' \code{trimming_threshold} (with smoothing via \code{smoothing_constant}) or applies overlap weighting 
#' if \code{overlap = TRUE}.  After estimating all propensity scores via cross‐fitting (in \code{cf_propensity}), 
#' it computes the density ratios and \(\phi\)-scores (in \code{compute_q_phi_score}), then proceeds
#' backward from $t = \tau$ to $t = 1$, at each step cross‐fitting a regression of the current pseudo‐outcome
#' on covariates to update the pseudo‐outcome for the next time step (in \code{cf_sequential_regression}).
#' The final pseudo‐outcome at time $t = 1$ serves as an influence‐function estimate for the overall flip effect.  
#'
#' If \code{life = TRUE}, the returned list also includes \code{$task} and \code{$control} so that 
#' one can call \code{flip_sdr_second_PO()} or \code{flip_sdr_treatment()} to continue estimation 
#' under a new target regime or to obtain a treatment‐level effect at a specific final time point, 
#' but using the same fold and propensity score estimate information.
#'
#' @examples
#' \dontrun{
#' # Basic usage with smooth trimming:
#' df <- sim_t4
#' df$A_1 <- rbinom(n = nrow(df), size = 1, prob = 0.3)
#' df$A_2 <- rbinom(n = nrow(df), size = 1, prob = 0.5)
#' df$A_3 <- rbinom(n = nrow(df), size = 1, prob = 0.7)
#' df$A_4 <- rbinom(n = nrow(df), size = 1, prob = 0.9)
#' df$X <- runif(n = nrow(df))
#'
#'flip_sdr(
#' data = df,
#' trt = c("A_1", "A_2", "A_3", "A_4"),
#' outcome = "Y",
#' baseline = c("X"),
#' time_vary = list(c("L_1"), c("L_2"), c("L_3"), c("L_4")),
#' cens               = NULL,
#' compete            = NULL,
#' id = "ID",
#' target_regime      = rep(1, 4),        # always treat
#' overlap            = FALSE,
#' trimming_threshold = 0.1,
#' smoothing_constant = 5,
#' outcome_type       = "continuous",
#' learners_outcome   = "SL.glm",
#' learners_trt       = "SL.glm",
#' folds              = 10
#' )
#' }
#'
#' @export
flip_sdr <- function(
    data,
    trt,
    outcome,
    baseline      = NULL,
    time_vary     = NULL,
    cens          = NULL,
    compete       = NULL,
    target_regime = NULL,
    overlap       = FALSE,
    trimming_threshold,
    smoothing_constant,
    k             = Inf,
    outcome_type  = c("binomial","continuous","survival"),
    id            = NULL,
    bounds        = NULL,
    learners_outcome = "SL.glm",
    learners_trt  = "SL.glm",
    folds         = 10,
    weights       = NULL,
    control       = lmtp_control(),
    life = FALSE
) {
  
  # Construct an lmtp task
  task <- LmtpTask$new(
    data = data,
    shifted = data,
    A = trt,
    Y = outcome,
    L = time_vary,
    W = baseline,
    C = cens,
    D = compete,
    k = k, 
    id = id,
    outcome_type = match.arg(outcome_type),
    folds = folds, 
    weights = weights,
    bounds = bounds
  )
  
  # Create progress bar object
  pb <- progressr::progressor(task$tau * folds * 2)
  
  # --- Estimate propensity scores --- #
  props <- cf_propensity(task, learners_trt, control, pb)
  
  # --- Define flip weight and derivative --- # 
  flip_funcs <- define_flip_functions(overlap, trimming_threshold, smoothing_constant)
  
  # --- Compute q scores, ratios, and phi scores --- #
  scores <- compute_q_phi_scores(
    task = task,
    prop_scores = props$prop_scores,
    target_regime    = target_regime,
    flip_func        = flip_funcs[[1]],
    flip_func_deriv  = flip_funcs[[2]]
  )

  # --- Run sequential regressions with cross-fitting --- #

  # Set initial pseudo-outcome equal to actual outcome
  task$natural$pseudo <- task$natural$final <- task$natural[, task$vars$Y]
  estims <- cf_sequential_regression(task, task$tau, scores, learners_outcome, control, pb)
  
  # Calculate final EIF
  ifvalues <- eif_Q(task, task$natural, estims$seq_reg_ests,
                    scores$q_scores, scores$phi_scores, scores$ratios, time = 1,
                    final_time = task$tau) 
    
    
  # Return outputs!
  out <- list(
    estimator = "SDR",
    ifvalues = ifvalues,
    estimate = ife::ife(x = mean(ifvalues),
                        eif = ifvalues - mean(ifvalues),
                        id = as.character(task$id)),
    prop_scores = props$prop_scores,
    fits_prop = props$fits,
    fits_m = estims$fits,
    outcome_type = task$outcome_type
  )
  
  if (life) {
    # Add task object to the output
    out$task <- task
    out$control <- control
  }
  
  class(out) <- "lmtp"
  
  return(out)
}

#' Internal: second‐stage SDR for potential outcomes
flip_sdr_second_PO <- function(flip,
                          target_regime,
                          overlap,
                          trimming_threshold,
                          smoothing_constant,
                          k,
                          outcome_type,
                          id,
                          bounds,
                          learners_outcome,
                          learners_trt,
                          folds,
                          weights,
                          control) {
  
  # Create progress bar object
  pb <- progressr::progressor(flip$task$tau * folds)
  
  # --- Define flip weight and derivative --- # 
  flip_funcs <- define_flip_functions(overlap, trimming_threshold, smoothing_constant)
  
  # --- Compute q scores, ratios, and phi scores --- #
  scores <- compute_q_phi_scores(
    task = flip$task,
    prop_scores = flip$prop_scores,
    target_regime    = target_regime,
    flip_func        = flip_funcs[[1]],
    flip_func_deriv  = flip_funcs[[2]]
  )
  
  # --- Run sequential regressions with cross-fitting --- #
  
  # Set initial pseudo-outcome equal to actual outcome
  flip$task$natural$pseudo <- flip$task$natural$final <- flip$task$natural[, flip$task$vars$Y]
  estims <- cf_sequential_regression(flip$task, flip$task$tau, scores, learners_outcome, control, pb)
  
  # Calculate final EIF
  ifvalues <- eif_Q(flip$task, flip$task$natural, estims$seq_reg_ests,
                    scores$q_scores, scores$phi_scores, scores$ratios, time = 1,
                    final_time = flip$task$tau) 
  # Return outputs
  out <- list(
    estimator = "SDR",
    ifvalues = ifvalues,
    estimate = ife::ife(x = mean(ifvalues),
                        eif = ifvalues - mean(ifvalues),
                        id = as.character(flip$task$id)),
    fits_m = estims$fits
    )
  
  class(out) <- "lmtp"
  
  return(out)
}

#' Internal: SDR for estimating average number of treatments
flip_sdr_treatment <- function(flip,
                             target_regime,
                             overlap,
                             trimming_threshold,
                             smoothing_constant,
                             k,
                             outcome_type,
                             id,
                             bounds,
                             learners_outcome,
                             folds,
                             weights,
                             control,
                             final_time) {
  
  # Progress bar
  pb <- progressr::progressor((final_time-1) * folds)
  
  # --- Define flip weight and derivative --- # 
  flip_funcs <- define_flip_functions(overlap, trimming_threshold, smoothing_constant)
  
  # --- Compute q scores, ratios, and phi scores --- #
  scores <- compute_q_phi_scores(
    task = flip$task,
    prop_scores = flip$prop_scores,
    target_regime    = target_regime,
    flip_func        = flip_funcs[[1]],
    flip_func_deriv  = flip_funcs[[2]]
  )
  
  # --- Run sequential regressions with cross-fitting --- #
  
  # Set initial pseudo-outcome equal to actual outcome
  flip$task$natural$pseudo <- flip$task$natural$final <- 
    pmax(0, pmin(1, scores$q_scores[, final_time] + scores$phi_scores[, final_time, 2]))
  
  if (final_time > 1) {
    estims <- cf_sequential_regression(flip$task, final_time-1, scores, learners_outcome, control, pb)
    fits_m <- estims$fits
    ifvalues <- eif_Q(flip$task, flip$task$natural, estims$seq_reg_ests,
                      scores$q_scores, scores$phi_scores, scores$ratios, time = 1,
                      final_time = final_time-1) 
  } else {
    fits_m <- vector("list", length = 0)
    ifvalues <- flip$task$natural$pseudo
  }

  
  # Return outputs
  out <- list(
    estimator = "SDR",
    ifvalues = ifvalues,
    estimate = ife::ife(x = mean(ifvalues),
                        eif = ifvalues - mean(ifvalues),
                        id = as.character(flip$task$id)),
    fits_m = fits_m
  )
  
  class(out) <- "lmtp"
  
  return(out)
}



