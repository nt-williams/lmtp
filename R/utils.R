get_folded_data <- function(x, ...) {
  UseMethod("get_folded_data")
}

#' @export
get_folded_data.data.frame <- function(x, folds, index, ...) {
  out <- list()
  out[["train"]] <- x[folds[[index]]$training_set, , drop = FALSE]
  out[["valid"]] <- x[folds[[index]]$validation_set, , drop = FALSE]
  out
}

#' @export
get_folded_data.matrix <- get_folded_data.data.frame

#' @export
get_folded_data.numeric <- function(x, folds, index, ...) {
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
    followed[, i] <- equal_elementwise(natural[[a]], shifted[[a]])
  }

  as.numeric(rowSums(!followed) == 0)
}

# Vectorised version of
#   mapply(function(x, y) isTRUE(all.equal(x, y)), as.list(x), as.list(y))
# Comparing one observation at a time with all.equal() makes followed_rule()
# the bottleneck of the density ratio estimation for large data. For numeric,
# factor and character vectors this returns the same result as the element-wise
# all.equal():
#   numeric:   relative difference <= tolerance (absolute when |x| <= tolerance)
#   factor:    same class and levels, then same value
#   character: same value
#   NA equals NA, and NA never equals a value.
# Other types fall back to the element-wise all.equal().
equal_elementwise <- function(x, y, tolerance = sqrt(.Machine$double.eps)) {
  if (is.numeric(x) && is.numeric(y)) {
    scale <- ifelse(is.finite(x) & abs(x) > tolerance, abs(x), 1)
    eq <- (x == y) | (abs(x - y) / scale <= tolerance)
  } else if (is.factor(x) && is.factor(y)) {
    if (!identical(class(x), class(y)) || !identical(levels(x), levels(y))) {
      return(rep(FALSE, length(x)))
    }
    eq <- as.character(x) == as.character(y)
  } else if (is.character(x) && is.character(y)) {
    eq <- x == y
  } else {
    return(mapply(function(x, y) isTRUE(all.equal(x, y)), as.list(x), as.list(y),
                  USE.NAMES = FALSE))
  }

  na_x <- is.na(x)
  na_y <- is.na(y)
  eq[na_x & na_y] <- TRUE
  eq[xor(na_x, na_y)] <- FALSE
  eq
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
  na_rows <- is.na(data[[x]])
  ohe <- model.matrix(~ -1 + as.character(data[[x]]))
  colnames(ohe) <- gsub("as.character\\(data\\[\\[x\\]\\]\\)", "", colnames(ohe))
  if (any(na_rows)) {
    full <- matrix(NA_real_, nrow = nrow(data), ncol = ncol(ohe), dimnames = list(NULL, colnames(ohe)))
    full[!na_rows, ] <- ohe
    ohe <- full
  }
  ohe
}

#' @export
summary.SuperLearner <- function(object, time = NULL, fold = NULL, level = NULL, ...) {
  values <- data.frame(risk = round(object$cvRisk, 3), coef = round(object$coef, 3))
  values <- data.frame(learner = rownames(values), values, row.names = NULL)
  data.table(time = time, fold = fold, level = level, values)
}

`%*0%` <- function(x, y) {
  res <- x * y
  res[is.na(res)] <- 0
  res
}