predict_delay_augment <- function(data, object, time, time_horizon,
                                  treatment, outcome, i, y1, d0, shifted = TRUE) {
  if (time > 1) {
    vars <- c(paste0("..i..lmtp_tmp_s", time - 1), treatment)
  } else {
    vars <- treatment
  }

  pred <- rep(NA_real_, nrow(data))

  if (!shifted) {
    pred[i] <- predict(object, data[i, ], 1e-05)
    return(calibrate(pred, y1, d0))
  }

  tmp <- data
  tmp[[treatment]] <- one_time_delay(data, vars)
  # If not the last time point, we replace s_t with a_t
  if (time < time_horizon) {
    tmp[[paste0("..i..lmtp_tmp_s", time)]] <- data[[treatment]]
  }

  pred[i] <- predict(object, tmp[i, ], 1e-05)
  calibrate(pred, y1, d0)
}

delay_augment <- function(x, ...) {
  UseMethod("delay_augment")
}

#' @export
delay_augment.data.frame <- function(x, sequences) {
  n <- nrow(sequences)
  t <- ncol(sequences)
  nr <- nrow(x)

  # Replicate data n times
  augmented <- x[rep(seq_len(nr), n), ]

  # Add all sequence columns at once
  seq_cols <- sequences[rep(seq_len(n), each = nr), , drop = FALSE]
  setnames(seq_cols, paste0("..i..lmtp_tmp_s", seq_len(t)))

  cbind(augmented, seq_cols)
}

#' @export
delay_augment.logical <- function(x, sequences) {
  n <- nrow(sequences)
  rep(x, n)
}

subset_augmented <- function(x, ...) {
  UseMethod("subset_augmented")
}

#' @export
subset_augmented.data.frame <- function(x, time, time_horizon) {
  if (time == time_horizon) return(x)
  id <- "..i..lmtp_id"
  if (time > 0) {
    id <- c(id, seqvars(seq_len(time)))
  }
  x[!collapse::fduplicated(x[, id, drop = FALSE]), ]
}

#' @export
subset_augmented.numeric <- function(x, id) {
  x[!collapse::fduplicated(id)]
}

# TODO: Need to figure out how to incorporate the censoring indicators
delay_riesz_rep <- function(data, treatments, propensity_scores, prob_observed,
                            this_treatment, this_censoring, this_time, time_horizon) {
  i <- seq_len(dim(propensity_scores)[1])

  # Create indicators (D in paper)
  ind <- rep(TRUE, nrow(data))
  for (time in this_time:time_horizon) {
    ind <- ind * one_time_delay(data, seqvars(this_time:time)) == data[[this_treatment]]
  }

  # Pre-allocate density ratios
  density_ratios <- rep(1, nrow(data))
  for (time in this_time:time_horizon) {
    column_obs <- match(as.character(data[[treatments[time]]]),
                        colnames(propensity_scores[, , time]))
    column_s <- match(as.character(data[[seqvars(time)]]),
                      colnames(propensity_scores[, , time]))
    prob_trt_natural <- propensity_scores[, , time][cbind(i, column_obs)]
    prob_trt_shifted <- propensity_scores[, , time][cbind(i, column_s)]
    density_ratios <- density_ratios * (prob_trt_shifted / prob_trt_natural)
  }

  ind %*0% density_ratios
}

delay_sdr_transformation <- function(data, pred_shifted, pred_natural,
                                     propensity_scores, prob_observed,
                                     trtvars, this_treatment, this_censoring,
                                     time, time_horizon, outcomevar, aug_key) {

  idvar <- "..i..lmtp_id"

  # Pre-compute the subset indices once (used for all k)
  subset_cols <- c(idvar, seqvars(seq_len(time - 1)))
  subset_keys <- aug_key[, subset_cols, drop = FALSE]

  # Vectorized accumulation over k = time:time_horizon
  pseudo <- 0
  for (k in time:time_horizon) {
    riesz <- delay_riesz_rep(data, trtvars, propensity_scores, prob_observed,
                             this_treatment, this_censoring, time, k)
    # Compute residual, take product with Riesz representers, then subset
    resid <- pred_shifted[[k + 1]] - pred_natural[[k]]
    pseudo <- pseudo + subset_augmented(riesz * resid, subset_keys)
  }

  pseudo
}
