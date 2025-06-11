#' @title Longitudinal Interventional Flip Effect SDR Estimator
#' @description
#' Estimate the longitudinal interventional flip effect comparing always‐treated versus never‐treated
#' using a sequentially doubly robust (SDR) procedure with either smooth trimming or overlap weighting.
#'
#' @param data               A \code{data.frame} containing all variables.
#' @param trt                Character vector of treatment column names (time‐varying, each ending in \code{"_tXX"}).
#' @param outcome            Name of the outcome column (character).
#' @param baseline           Character vector of baseline covariate names, or \code{NULL}.
#' @param time_vary          Character vector of time‐varying covariate names (each ending in \code{"_tXX"}), or \code{NULL}.
#' @param cens               Name of the censoring indicator column (character), or \code{NULL}.
#' @param compete            Name of the competing event indicator column (character), or \code{NULL}.
#' @param overlap            Logical; if \code{TRUE}, use overlap weights at each time point, otherwise use smooth trimming.
#' @param trimming_threshold Numeric threshold for trimming (or weight cutoff) when \code{overlap = FALSE}.
#' @param smoothing_constant Numeric smoothing parameter (larger → weight function closer to a hard indicator).
#' @param k                  Integer (or \code{Inf}); history‐window size for the underlying \code{LmtpTask}.
#'                           If \code{Inf}, use the full treatment history.
#' @param outcome_type       Character; one of \code{c("binomial","continuous","survival")}.
#'                           Determines whether pseudo‐outcomes are bounded.
#' @param id                 Name of the subject‐identifier column (character), or \code{NULL}.
#' @param bounds             Optional numeric vector of length 2 (\code{c(lower, upper)}) giving outcome bounds.
#'                           If \code{NULL}, no bounding is applied on any regression fits.
#' @param learners_outcome   Character vector of SuperLearner (or other) learners for sequential outcome regressions.
#' @param learners_trt       Character vector of SuperLearner (or other) learners for time‐varying treatment models.
#' @param folds              Integer ≥ 2; number of cross‐fitting folds for propensity and outcome regressions.
#' @param weights            Optional numeric vector of sample weights (length = \code{nrow(data)}).
#'                           If \code{NULL}, no sample weighting is applied.
#' @param control            An \code{lmtp_control} object controlling learner options and parallelization.
#' @param num_timepoints     Integer; total number of time points (i.e.\ \(\tau\)). Used for loops over times.
#'
#' @return A list with class \code{life}, containing:
#' \describe{
#'   \item{\code{estimator}}{Character; the string \code{"SDR"}, indicating the estimator used.}
#'   \item{\code{ifvalues}}{Numeric vector of length \(n\), giving the final influence‐function values
#'                          (ratio of always‐ vs never‐treated flip effects).}
#'   \item{\code{estimate}}{An \code{\link[ife]{ife}} object with \code{x =} point estimate of
#'                          the flip‐effect ratio and \code{eif = ifvalues}.}
#'   \item{\code{auxiliary_info}}{A list containing:
#'     \itemize{
#'       \item{\code{always_treated_avg_PO}}{Result of \code{flip_sdr} under always‐treated regime.}
#'       \item{\code{never_treated_avg_PO}}{Result of \code{flip_sdr_second_PO} under never‐treated regime.}
#'       \item{\code{always_treated_avg_tx}}{List of \code{flip_sdr_treatment} outputs for each time \(t\) under always‐treated.}
#'       \item{\code{never_treated_avg_tx}}{List of \code{flip_sdr_treatment} outputs for each time \(t\) under never‐treated.}
#'     }
#'   }
#' }
#'
#' @details
#' This function estimates the longitudinal interventional flip effects of 
#' (1) always‐treated versus (2) never‐treated. Internally:
#' \enumerate{
#'   \item It calls \code{flip_sdr(..., target_regime = rep(1, num_timepoints), life = TRUE)} 
#'         to estimate the average potential outcome when everyone is flipped towards treatment
#'         at all times (via smooth trimming or overlap weighting) and stores the resulting
#'         \code{flip} object as \code{always_treated_avg_PO}.
#'   \item It loops backward through \(\{\,\tau, \dots, 1\,\}\), calling
#'         \code{flip_sdr_treatment(always_treated_avg_PO, final_time = t)} to compute the
#'         average number of treatments at each time \(t\) under always‐treated; stored in
#'         \code{always_treated_avg_tx}.
#'   \item It calls \code{flip_sdr_second_PO(always_treated_avg_PO, target_regime = rep(0, num_timepoints))} 
#'         to estimate the average potential outcome when everyone is flipped towards control
#'         (smooth trimming or overlap weighting under never treated); stored as \code{never_treated_avg_PO}.
#'   \item It loops backward again, calling
#'         \code{flip_sdr_treatment(always_treated_avg_PO, target_regime = rep(0, num_timepoints), final_time = t)}
#'         to compute the average number of treatments at each time \(t\) under never‐treated;
#'         stored in \code{never_treated_avg_tx}.
#'   \item Finally, it returns influence function estimates for the longitudinal interventional
#'         flip effect alongside auxiliary info from previous steps, including the estimates 
#'         for the difference in number of treatments at each timepoint
#' }
#' 
#' @examples
#' \dontrun{
#' # Suppose you have 4 time points and want to compare always vs never treat:
#' mydata <- sim_t4
#' mydata$A_1 <- rbinom(n = nrow(mydata), size = 1, prob = 0.3)
#' mydata$A_2 <- rbinom(n = nrow(mydata), size = 1, prob = 0.5)
#' mydata$A_3 <- rbinom(n = nrow(mydata), size = 1, prob = 0.7)
#' mydata$A_4 <- rbinom(n = nrow(mydata), size = 1, prob = 0.9)
#' mydata$X <- runif(n = nrow(mydata))
#' life_sdr(
#'   data = mydata, 
#'   trt = c("A_1", "A_2", "A_3", "A_4"),
#'   outcome = "Y",
#'   baseline = c("X"),
#'   time_vary = list(c("L_1"), c("L_2"), c("L_3"), c("L_4")),
#'   cens               = NULL,
#'   compete            = NULL,
#'   id = "ID",
#'   overlap            = FALSE,
#'   trimming_threshold = 0.1,
#'   smoothing_constant = 5,
#'   outcome_type       = "continuous",
#'   learners_outcome   = "SL.glm",
#'   learners_trt       = "SL.glm",
#'   folds              = 10,
#'   num_timepoints = 4
#'   )
#' }
#'
#' @export
life_sdr <- function(
    data,
    trt,
    outcome,
    baseline      = NULL,
    time_vary     = NULL,
    cens          = NULL,
    compete       = NULL,
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
    num_timepoints
) {
  
  always_treated <- flip_sdr(
    data = data, 
    trt = trt, 
    outcome = outcome,
    baseline = baseline,
    time_vary = time_vary,
    cens = cens,
    compete = compete,
    id = id,
    target_regime = rep(1, num_timepoints), # always treat
    overlap = overlap,
    trimming_threshold = trimming_threshold,
    smoothing_constant = smoothing_constant,
    outcome_type       = outcome_type,
    learners_outcome   = learners_outcome,
    learners_trt       = learners_trt,
    folds              = folds,
    life = T
  )
  
  # Construct estimates for average number of treatments under always treated
  # Re-use information from first run (e.g., propensity score estimates, 
  # fold assignments, id, etc.)
  always_avg_treatments <- vector("list", length = num_timepoints) 
  for (time in num_timepoints:1) {
    always_avg_treatments[[time]] <-
      flip_sdr_treatment(always_treated,
                       target_regime = rep(1, num_timepoints),
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
                       final_time = time) 
  }
  
  # Estimate mean potential outcome (and re-use information from first run)
  never_treated <- flip_sdr_second_PO(
    flip = always_treated,
    target_regime = rep(0, num_timepoints), # never treated
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
    control) 
  
  # Construct estimates for average number of treatments under never treated
  #  (and re-use information from first run)
  never_avg_treatments <- vector("list", length = num_timepoints) 
  for (time in num_timepoints:1) {
    never_avg_treatments[[time]] <-
      flip_sdr_treatment(always_treated,
                       target_regime = rep(0, num_timepoints),
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
                       final_time = time) 
  }

  # Compute numerator IFs and point estimate (always_avg_PO - never_avg_PO)
  num_ifs <- always_treated$ifvalues - never_treated$ifvalues
  num_est <- mean(num_ifs)

  # Compute denominator: average per-timepoint absolute change in treatment
  trt_diff_info <- data.frame()
  denom_est <- 0
  denom_ifs <- rep(0, length(id))
  for (t in 1:always_treated$task$tau) {
    
    always_ifvalues <- always_avg_treatments[[t]]$ifvalues
    never_ifvalues  <- never_avg_treatments[[t]]$ifvalues
    
    this_time_est <- abs(mean(always_ifvalues - never_ifvalues)) 
    this_time_ifs <- sign(mean(always_ifvalues - never_ifvalues)) * 
      (always_ifvalues - never_ifvalues)
    
    denom_est <- denom_est + this_time_est / always_treated$task$tau
    denom_ifs <- denom_ifs + this_time_ifs / always_treated$task$tau
    
    trt_diff_info <- rbind(trt_diff_info, 
      data.frame(
        t = t,
        ptest = this_time_est,
        sdest = sd(this_time_ifs) / sqrt(length(this_time_ifs))
      )
    )
  }
  
  # Form ratio influence‐function and return
  ratio_ifs <- num_ifs / denom_est - (num_est / denom_est^2) * denom_ifs 
  
  return(
    list(
      estimator = "SDR",
      ifvalues = ratio_ifs,
      estimate =   ife::ife(x = num_est / denom_est,
                            eif = ratio_ifs,
                            id = as.character(always_treated$task$id)),
      auxiliary_info = list(
        trt_diff_info = trt_diff_info,
        always_treated_avg_PO = always_treated,
        never_treated_avg_PO = never_treated,
        always_treated_avg_tx = always_avg_treatments,
        never_treated_avg_tx = never_avg_treatments,
        overall_num_ifs = num_ifs,
        overall_denom_ifs = denom_ifs
      )
    )
  )
}


