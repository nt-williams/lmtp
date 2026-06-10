get_folded_data <- function(x, ...) {
  UseMethod("get_folded_data")
}

#' @export
get_folded_data.data.frame <- function(x, folds, index) {
  out <- list()
  out[["train"]] <- x[folds[[index]]$training_set, , drop = FALSE]
  out[["valid"]] <- x[folds[[index]]$validation_set, , drop = FALSE]
  out
}

#' @export
get_folded_data.matrix <- get_folded_data.data.frame

#' @export
get_folded_data.numeric <- function(x, folds, index) {
  out <- list()
  out[["train"]] <- x[folds[[index]]$training_set]
  out[["valid"]] <- x[folds[[index]]$validation_set]
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

  followed <- matrix(nrow = nrow(natural), ncol = length(A))

  for (i in seq_along(A)) {
    a <- A[i]
    followed[, i] <- mapply(function(x, y) isTRUE(all.equal(x, y)),
                            as.list(natural[, a]),
                            as.list(shifted[, a]))
  }

  apply(followed, 1, prod)
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

`%and%` <- function(o, r) {
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

current_trt <- function(trt, time) {
  if (length(trt) > 1) {
    return(trt[[time]])
  }
  trt[[1]]
}

calibrate <- function(pred, prior_free, comp_free) {
  pred[!prior_free] <- 0
  pred[!comp_free] <- 1
  pred
}

one_hot_encode <- function(data, x) {
  ohe <- model.matrix(~ -1 + as.character(data[[x]]))
  colnames(ohe) <- gsub("as.character\\(data\\[\\[x\\]\\]\\)", "", colnames(ohe))
  ohe
}

#' @export
summary.SuperLearner <- function(object, time = NULL, fold = NULL, level = NULL, ...) {
  values <- data.frame(risk = round(object$cvRisk, 3), coef = round(object$coef, 3))
  values <- data.frame(learner = rownames(values), values, row.names = NULL)
  data.table(time = time, fold = fold, level = level, values)
}
