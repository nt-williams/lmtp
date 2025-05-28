#' @keywords internal
#' @title Cross‐fitted propensity score estimation
#' @description
#' Perform V‐fold cross‐validation to estimate time‐varying propensity scores.
#'
#' @param task     An `lmtp_task` object created by `lmtp_task()`.
#' @param learners Character; either "rf" (ranger) or "sl" (SuperLearner).
#' @param control  An `lmtp_control()` object; supplies `.learners_trt_folds` and `.quiet`.
#' @param pb       Function; optional progress callback called after each timepoint.
#'
#' @return A list with:
#'   \item{prop_scores}{Array of dimension (n_valid × T × 2): P(A_t=0) in [,,1], P(A_t=1) in [,,2].}
#'   \item{fits}{List of length T with fitted model objects or weights.}
cf_propensity <- function(task, learners, control, pb = NULL) {
  folds <- task$folds
  nf    <- length(folds)
  ans   <- vector("list", length = nf)
  
  for (f in seq_len(nf)) {
    ans[[f]] <- future::future({
      estimate_propensity(task, f, learners, control, pb)
    }, seed = TRUE)
  }
  ans <- future::value(ans)
  
  prop_list <- lapply(ans, function(x) x$prop_scores)
  prop_scores <- recombine(rbind_depth(prop_list, NULL), folds)
  fits <- lapply(ans, function(x) x$fits)
  
  list(prop_scores = prop_scores, fits = fits)
}

#' @keywords internal
#' @title Fold‐specific propensity score estimation
#' @description
#' Estimate time‐varying propensity scores within a single cross‐fit fold, 
#' fitting only on at‐risk subjects.
#'
#' @inheritParams cf_propensity
#' @param fold Integer index of the fold to use as evaluation set.
#'
#' @return A list with:
#'   \item{prop_scores}{Array (n_valid × T × 2) of estimated propensities on validation set.}
#'   \item{fits}{List of length T of fitted models or weight vectors.}
estimate_propensity <- function(task, fold, learners, control, pb = NULL) {
  
  # Pull natural train/validation splits
  natural    <- get_folded_data(task$natural, task$folds, fold)
  train_data <- natural$train
  valid_data <- natural$valid
  
  n_valid <- nrow(valid_data)
  # Only estimate for timepoints with treatment data
  T_trt   <- length(task$vars$A)
  
  # Initialize outputs
  prop_scores <- array(NA_real_, dim = c(n_valid, T_trt, 2))
  fits        <- vector("list", T_trt)
  
  baseline <- task$baseline
  id_var   <- task$vars$id
  
  for (t in seq_len(T_trt)) {
    
    # Treatment at time t
    A_t <- task$vars$A[[t]]
    
    # At-risk subset indices
    train_idx <- ii(task$observed(train_data, t - 1), task$R(train_data, t))
    valid_idx <- ii(task$observed(valid_data, t - 1), task$R(valid_data, t))
    
    # Predictors: baseline + treatment history up to t-1 + time‐varying covariates up to t
    pred_cols <- c(
      baseline,
      task$vars$history("A", t),
      task$vars$history("L", t + 1)
    )
    
    # Training data frame
    train_df <- train_data[train_idx, c(id_var, pred_cols, A_t), drop = FALSE]
    
    # Fit propensity model
    fit <- run_ensemble(
      data     = train_df,
      outcome  = A_t,
      learners = learners,
      family   = "binomial",
      id       = id_var,
      folds    = control$.learners_trt_folds
    )
    fits[[t]] <- if (control$.return_full_fits) fit else extract_sl_weights(fit)
    
    # Validation data frame and predictions
    valid_df <- valid_data[valid_idx, c(id_var, pred_cols, A_t), drop = FALSE]
    preds    <- predict(fit, valid_df)
    
    # Store P(A_t=1 | H_t) and P(A_t=0 | H_t)
    prop_scores[valid_idx, t, 2] <- preds
    prop_scores[valid_idx, t, 1] <- 1 - preds
    
    if (!is.null(pb)) pb()
  }
  
  list(prop_scores = prop_scores, fits = fits)
}