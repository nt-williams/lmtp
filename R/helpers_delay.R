predict_delay_augment <-
  function(data,
           object,
           sequences,
           time,
           time_horizon,
           treatment,
           outcome,
           i, y1, d0,
           shifted = TRUE) {
    if (time > 1) {
      vars <- c(paste0("..i..lmtp_tmp_s", time - 1), treatment)
    } else {
      vars <- treatment
    }

    tmp <- data
    tmp[[treatment]] <- one_time_delay(data, vars)
    # If not the last time point, we replace s_t with a_t
    if (time < time_horizon) {
      tmp[[paste0("..i..lmtp_tmp_s", time)]] <- data[[treatment]]
    }
    ai <- delay_augment(i, sequences)

    data[[outcome]] <- NA_real_
    data[ai, outcome] <- predict(object, tmp[ai, ], 1e-05)

    # Handle deterministic predictions from survival/competing risks
    data[which(!delay_augment(y1, sequences)), outcome] <- 0
    data[which(!delay_augment(d0, sequences)), outcome] <- 1

    data
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
