check_lmtp_data <- function(x, trt, outcome, baseline, time_vary, cens, id) {
  for (t in 1:determine_tau(outcome, trt)) {
    ci <- censored(x, cens, t)$j
    di <- at_risk(x, risk_indicators(outcome), t, TRUE)
    if (length(trt) > 1) {
      trt_t <- trt[[t]]
    } else {
      trt_t <- trt[[1]]
    }
    data_t <- x[ci & di, c(trt_t, baseline, unlist(time_vary[t])), drop = FALSE]

    if (any(is.na(data_t))) {
      return("Missing data found in treatment and/or covariate nodes for uncensored observations")
    }
  }

  TRUE
}

assertLmtpData <- checkmate::makeAssertionFunction(check_lmtp_data)

assert_trt <- function(trt, tau) {
  is_list <- is.list(trt)
  if (!isTRUE(is_list)) {
    return(assertTrtCharacter(trt, tau))
  }
  assertTrtList(trt, tau)
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

assertTrtCharacter <- checkmate::makeAssertionFunction(check_trt_character)

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

assertTrtList <- checkmate::makeAssertionFunction(check_trt_list)

check_reserved_names <- function(x) {
  bad_names <- c("lmtp_id", "tmp_lmtp_stack_indicator", "tmp_lmtp_scaled_outcome") %in% names(x)
  if (!any(bad_names)) {
    return(TRUE)
  }
  "'lmtp_id', 'tmp_lmtp_stack_indicator', and 'tmp_lmtp_scaled_outcome' are reserved variable names"
}

assertReservedNames <- checkmate::makeAssertionFunction(check_reserved_names)

check_shifted_data <- function(x, natural, doesnt_change, cens, null.ok = TRUE) {
  if (is.null(x)) {
    if (null.ok)
      return(TRUE)
    return("Can't be 'NULL'")
  }

  if (!(identical(natural[doesnt_change], x[doesnt_change]))) {
    return("The only columns that can be different between `data` and `shifted` are those indicated in `trt` and `cens`")
  }

  if (is.null(cens)) {
    return(TRUE)
  }

  if (!all(x[cens] == 1)) {
    return("Censoring variables should be 1 in 'shifted'")
  }

  TRUE
}

assertShiftedData <- checkmate::makeAssertionFunction(check_shifted_data)

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

assert_not_data_table <- assertNotDataTable <- checkmate::makeAssertionFunction(check_not_data_table)

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

assertOutcomeTypes <- checkmate::makeAssertionFunction(check_outcome_types)

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

assertContrastType <- checkmate::makeAssertionFunction(check_contrast_type)

check_lmtp_list <- function(x) {
  all_lmtp <- all(unlist(lapply(x, is.lmtp)))
  if (!all_lmtp) {
    return("Objects must be of type 'lmtp'")
  }
  TRUE
}

assertLmtpList <- checkmate::makeAssertionFunction(check_lmtp_list)

check_dr <- function(x) {
  all_dr <- all(lapply(x, function(x) x[["estimator"]]) %in% c("TMLE", "SDR"))
  if (!all_dr) {
    return("Contrasts not implemented for substitution/IPW estimators")
  }
  TRUE
}

assertDr <- checkmate::makeAssertionFunction(check_dr)

check_ref_class <- function(x) {
  if (!is.lmtp(x)) {
    is_num <- checkmate::test_number(x)
    if (isFALSE(is_num)) {
      return("Must either be a single numeric value or another lmtp object")
    }
  }
  TRUE
}

assertRefClass <- checkmate::makeAssertionFunction(check_ref_class)

check_trt_type <- function(data, trt, mtp) {
  is_decimal <- vector("logical", length(trt))
  for (i in seq_along(trt)) {
    a <- data[[trt[i]]]
    if (is.character(a) | is.factor(a)) next
    is_decimal[i] <- any(schoolmath::is.decimal(a[!is.na(a)]))
  }
  if (any(is_decimal) & isFALSE(mtp)) {
      cli::cli_warn("Detected decimalish `trt` values and {.code mtp = FALSE}. Consider setting {.code mtp = TRUE} if getting errors.")
  }
}

check_same_weights <- function(weights) {
  if (length(weights) == 1) {
    check <- TRUE
  } else if (length(weights) == 2) {
    check <- identical(weights[[1]], weights[[2]])
  } else {
    check <- all(sapply(1:(length(weights) - 1), function(i) identical(weights[[i]], weights[[i + 1]])))
  }

  if (isFALSE(check)) {
    return("Weights must all be the same.")
  }
  TRUE
}

assertSameWeights <- checkmate::makeAssertionFunction(check_same_weights)

