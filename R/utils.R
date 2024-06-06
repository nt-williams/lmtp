determine_tau <- function(outcome, trt) {
  surv <- length(outcome) > 1
  if (!surv) {
    return(length(trt))
  }
  length(outcome)
}

setup_cv <- function(data, V = 10, id, strata, outcome_type) {
  if (length(unique(id)) == nrow(data) & outcome_type == "binomial") {
    strata <- data[[strata]]
    strata[is.na(strata)] <- 2
    out <- origami::make_folds(data, V = V, strata_ids = strata)
  } else {
    out <- origami::make_folds(data, cluster_ids = id, V = V)
  }

  if (V > 1) {
    return(out)
  }
  out[[1]]$training_set <- out[[1]]$validation_set
  out
}

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

scale_y <- function(y, bounds) {
  if (is.null(bounds)) {
    return(y)
  }
  (y - bounds[1]) / (bounds[2] - bounds[1])
}

y_bounds <- function(y, outcome_type, bounds = NULL) {
  if (outcome_type == "binomial" || is.null(outcome_type)) {
    return(NULL)
  }
  if (is.null(bounds)) {
    return(c(min(y, na.rm = T), max(y, na.rm = T)))
  }
  c(bounds[1], bounds[2])
}

rescale_y_continuous <- function(scaled, bounds) {
  (scaled*(bounds[2] - bounds[1])) + bounds[1]
}

censored <- function(data, cens, tau) {
  # when no censoring return TRUE for all obs
  if (is.null(cens)) {
    return(list(i = rep(TRUE, nrow(data)), j = rep(TRUE, nrow(data))))
  }

  # other wise find censored observations
  i <- data[[cens[tau]]] == 1

  if (tau > 1) {
    return(list(i = i, j = data[[cens[tau - 1]]] == 1))
  }

  list(i = i, j = rep(TRUE, nrow(data)))
}

at_risk <- function(data, risk, tau, check = FALSE) {
  if (is.null(risk)) {
    return(rep(TRUE, nrow(data)))
  }

  if (tau == 1) {
    return(rep(TRUE, nrow(data)))
  }

  if (check) {
    return(data[[risk[tau - 1]]] == 0 & !is.na(data[[risk[tau - 1]]]))
  }

  data[[risk[tau - 1]]] == 1 & !is.na(data[[risk[tau - 1]]])
}

followed_rule <- function(obs_trt, shifted_trt, mtp) {
  if (mtp) {
    if (inherits(obs_trt, "data.frame")) {
      return(rep(TRUE, nrow(obs_trt)))
    }
    return(rep(TRUE, length(obs_trt)))
  }

  mapply(function(x, y) isTRUE(all.equal(x, y)), as.list(obs_trt), as.list(shifted_trt))
}

transform_sdr <- function(r, tau, max, shifted, natural) {
  natural[is.na(natural)] <- -999
  shifted[is.na(shifted)] <- -999
  m <- shifted[, (tau + 2):(max + 1), drop = FALSE] - natural[, (tau + 1):max, drop = FALSE]
  rowSums(r * m, na.rm = TRUE) + shifted[, tau + 1]
}

recombine_ratios <- function(x, folds) {
  ind <- Reduce(c, lapply(folds, function(x) x[["validation_set"]]))

  returns <- list()

  returns$ratios <- Reduce(
    rbind,
    lapply(x, function(x) x[["ratios"]])
  )[order(ind), ]

  if (is.null(dim(returns[["ratios"]]))) {
    returns[["ratios"]] <- as.matrix(
      returns[["ratios"]],
      nrow = length(returns[["ratios"]]),
      ncol = 1
    )
  }

  returns$fits <- lapply(x, function(x) x[["fits"]])
  returns
}

trim_ratios <- function(x, trim) {
  x[["ratios"]] <- pmin(x[["ratios"]], quantile(x[["ratios"]], trim))
  x
}

recombine_outcome <- function(x, part, folds) {
  ind <- Reduce(c, lapply(folds, function(x) x[["validation_set"]]))
  Reduce(rbind, lapply(x, function(x) x[[part]]))[order(ind), , drop = FALSE]
}

is.lmtp <- function(x) {
  class(x) == "lmtp"
}

sw <- function(x) {
  suppressWarnings(x)
}

final_outcome <- function(outcomes) {
  outcomes[length(outcomes)]
}

extract_sl_weights <- function(fit) {
  fit$coef
}

#' Time To Event Last Outcome Carried Forward
#'
#' A helper function to prepare survival data for use with LMTP estimators
#' by imputing outcome nodes using last outcome carried forward when an observation
#' experiences the event before the end-of-follow-up.
#'
#' @param data The dataset to modify.
#' @param outcomes A vector of outcome nodes ordered by time.
#'
#' @return A modified dataset with future outcome nodes set to 1 if an observation
#'  experienced an event at any previous time point.
#'
#' @importFrom data.table as.data.table `:=` .SD
#' @export
#' @examples
#' event_locf(sim_point_surv, paste0("Y.", 1:6))
event_locf <- function(data, outcomes) {
  DT <- as.data.table(data)
  tau <- length(outcomes)
  for (j in outcomes[1:(tau - 1)]) {
    modify <- setdiff(outcomes[match(j, outcomes):tau], j)
    DT[get(j) == 1 & !is.na(get(j)), (modify) := lapply(.SD, function(x) 1), .SDcols = modify]
  }
  DT[]
  DT
}

create_ids <- function(data, id) {
  if (is.null(id)) {
    return(1:nrow(data))
  }
  data[[id]]
}

convert_to_surv <- function(x) {
  data.table::fcase(
    x == 0, 1,
    x == 1, 0
  )
}

missing_outcome <- function(x) {
  ifelse(is.na(x), 0, x)
}

risk_indicators <- function(x) {
  if (length(x) == 1) {
    return(NULL)
  }

  x[1:(length(x) - 1)]
}

compute_weights <- function(r, t, tau) {
  out <- t(apply(r[, t:tau, drop = FALSE], 1, cumprod))
  if (ncol(out) > ncol(r)) return(t(out))
  out
}

is_normalized <- function(x, tolerance = .Machine$double.eps^0.5) {
  # Check if the mean is approximately 1 within the given tolerance
  abs(mean(x) - 1) < tolerance
}
