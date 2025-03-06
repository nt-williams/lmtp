get_folded_data <- function(data, folds, index) {
  out <- list()
  out[["train"]] <- data[folds[[index]]$training_set, , drop = FALSE]
  out[["valid"]] <- data[folds[[index]]$validation_set, , drop = FALSE]
  out
}

fix_censoring_ind <- function(data, cens) {
  if (is.null(cens)) {
    return(data)
  }

  data <- data.table::copy(data)
  for (cen in cens) {
    data.table::set(data, j = cen, value = ifelse(is.na(data[[cen]]), 0, data[[cen]]))
  }
  data
}

bound <- function(x, p = 1e-05) {
  pmax(pmin(x, 1 - p), p)
}

followed_rule <- function(natural, shifted, A, mtp) {
  if (mtp) {
    return(rep(TRUE, nrow(natural)))
  }

  mapply(function(x, y) isTRUE(all.equal(x, y)), as.list(natural[, A]), as.list(shifted[, A]))
}

trim <- function(x, trim) {
  pmin(x, quantile(x, trim, na.rm = TRUE))
}

is.lmtp <- function(x) {
  class(x) == "lmtp"
}

sw <- function(x) {
  suppressWarnings(x)
}

last <- function(x) {
  x[length(x)]
}

extract_sl_weights <- function(fit) {
  if (inherits(fit, "mlr3superlearner")) {
    return(cbind(Risk = fit$risk))
  }
  fit$coef
}

convert_to_surv <- function(x) {
  data.table::fcase(
    x == 0, 1,
    x == 1, 0
  )
}

is_normalized <- function(x, tolerance = .Machine$double.eps^0.5) {
  # Check if the mean is approximately 1 within the given tolerance
  abs(mean(x) - 1) < tolerance
}

fix_surv_time1 <- function(x) {
  to_fix <- x[[1]]$estimate
  x[[1]]$estimate <- ife::ife(1 - to_fix@x, 1 - to_fix@eif, to_fix@weights, to_fix@id)
  x
}

is_decimal <- function(x) {
  test <- floor(x)
  !(x == test)
}

ii <- function(o, r) {
  i <- vector("logical", length(o))
  for (j in 1:length(o)) {
    if (is.na(r[j]) & !is.na(o[j])) {
      i[j] <- o[j]
    } else if (!is.na(r[j]) & is.na(o[j])) {
      i[j] <- r[j]
    } else {
      i[j] <- o[j] & r[j]
    }
  }
  i
}
