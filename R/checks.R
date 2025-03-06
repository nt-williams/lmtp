check_lmtp_data = function(self) {
  for (t in 1:self$tau) {
    r <- self$R(self$natural, t)
    o <- self$observed(self$natural, t - 1)
    i <- ii(o, r)

    if (length(self$vars$A) > 1) {
      A_t <- self$vars$A[[t]]
    } else {
      A_t <- self$vars$A[[1]]
    }

    data_t <- self$natural[i, c(A_t, self$vars$W, unlist(self$vars$L[t])), drop = FALSE]

    if (any(is.na(data_t))) {
      return("Missing data found in treatment and/or covariate nodes for uncensored observations")
    }
  }

  TRUE
}

assert_lmtp_data <- checkmate::makeAssertionFunction(check_lmtp_data)

assert_trt <- function(trt, tau) {
  is_list <- is.list(trt)
  if (!isTRUE(is_list)) {
    return(assert_trt_character(trt, tau))
  }
  assert_trt_list(trt, tau)
}

check_trt_character <- function(trt, tau) {
  is_character <- checkmate::check_character(trt)
  if (!isTRUE(is_character)) {
    return(is_character)
  }

  if (length(trt) != 1 && length(trt) != tau) {
    return(paste0("'trt' should be of length 1 or ", tau))
  }

  TRUE
}

assert_trt_character <- checkmate::makeAssertionFunction(check_trt_character)

check_trt_list <- function(trt, tau) {
  is_list <- checkmate::check_list(trt)
  if (!isTRUE(is_list)) {
    return(is_list)
  }

  if (length(trt) != 1 && length(trt) != tau) {
    return(paste0("'trt' should be of length 1 or ", tau))
  }

  TRUE
}

assert_trt_list <- checkmate::makeAssertionFunction(check_trt_list)

check_shifted_data <- function(natural, shifted, trt, cens) {
  is_same <- setdiff(names(natural), c(trt, cens))

  if (!(identical(natural[is_same], shifted[is_same]))) {
    return("The only columns that can be different between `data` and `shifted` are those indicated in `trt` and `cens`")
  }

  if (is.null(cens)) {
    return(TRUE)
  }

  if (!all(shifted[cens] == 1)) {
    return("Censoring variables should be 1 in 'shifted'")
  }

  TRUE
}

assert_correctly_shifted <- checkmate::makeAssertionFunction(check_shifted_data)

check_not_data_table <- function(x) {
  is_data_frame <- checkmate::checkDataFrame(x)
  if (!isTRUE(is_data_frame)) {
    return(is_data_frame)
  }

  is_data_table <- data.table::is.data.table(x)
  if (is_data_table) {
    return("Must be a 'data.frame', not a 'data.table'")
  }
  TRUE
}

assert_not_data_table <- checkmate::makeAssertionFunction(check_not_data_table)

check_outcome_types <- function(x, outcomes, outcome_type) {
  x <- x[, outcomes, drop = FALSE]
  all_numeric <- checkmate::testDataFrame(x, types = "numeric")
  if (!all_numeric) {
    return("Outcome variables must be of type numeric")
  }

  if (outcome_type %in% c("binomial", "survival")) {
    vals <- lapply(outcomes, function(var) as.character(unique(na.omit(x[[var]]))))
    all_binary <- all(unlist(vals) %in% c("0", "1"))

    if (!isTRUE(all_binary))
      return("Only 0 and 1 allowed in outcome variables if 'outcome_type' set to binomial or survival")
  }
  TRUE
}

assert_outcome_types <- checkmate::makeAssertionFunction(check_outcome_types)

check_contrast_type <- function(x, fits) {
  if (x == "additive") {
    return(TRUE)
  }

  all_binom <- all(lapply(fits, function(x) x[["outcome_type"]]) == "binomial")
  if (!all_binom) {
    return(paste0("'", x, "' specified but one or more outcome types are not 'binomial' or 'survival'"))
  }

  TRUE
}

assert_contrast_type <- checkmate::makeAssertionFunction(check_contrast_type)

check_lmtp_list <- function(x) {
  all_lmtp <- all(unlist(lapply(x, is.lmtp)))
  if (!all_lmtp) {
    return("Objects must be of type 'lmtp'")
  }
  TRUE
}

assert_lmtp_list <- checkmate::makeAssertionFunction(check_lmtp_list)

check_ref_class <- function(x) {
  if (!is.lmtp(x)) {
    is_num <- checkmate::test_number(x)
    if (isFALSE(is_num)) {
      return("Must either be a single numeric value or another lmtp object")
    }
  }
  TRUE
}

assert_ref_class <- checkmate::makeAssertionFunction(check_ref_class)

check_trt_type <- function(data, trt, mtp) {
  is_decimal <- vector("logical", length(trt))
  for (i in seq_along(trt)) {
    a <- data[[trt[i]]]
    if (is.character(a) | is.factor(a)) next
    is_decimal[i] <- any(is_decimal(a[!is.na(a)]))
  }
  if (any(is_decimal) & isFALSE(mtp)) {
      cli::cli_warn("Detected decimalish `trt` values and {.code mtp = FALSE}. Consider setting {.code mtp = TRUE} if getting errors.")
  }
}
