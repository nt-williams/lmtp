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

delay_augment.logical <- function(x, sequences) {
  n <- nrow(sequences)
  rep(x, n)
}

subset_augmented <- function(x, ...) {
  UseMethod("subset_augmented")
}

subset_augmented.data.frame <- function(x, time, time_horizon) {
  if (time == time_horizon) return(x)
  id <- "..i..lmtp_id"
  if (time > 0) {
    id <- c(id, seqvars(seq_len(time)))
  }
  x[!collapse::fduplicated(x[, id, drop = FALSE]), ]
}

subset_augmented.numeric <- function(x, id) {
  x[!collapse::fduplicated(id)]
}

delay_riesz_rep <- function(data, treatments, propensity_scores,
                            this_treatment, this_time, time_horizon) {
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
