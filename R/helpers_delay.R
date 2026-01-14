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
