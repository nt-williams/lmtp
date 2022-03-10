check_lmtp_data <- function(x, trt, outcome, baseline, time_vary, cens, id) {
  is_data_frame <- checkmate::checkDataFrame(x)
  if (!isTRUE(is_data_frame)) {
    return(is_data_frame)
  }

  good_outcome <- checkmate::checkCharacter(outcome)
  if (!isTRUE(good_outcome)) {
    return(good_outcome)
  }

  good_trt <- checkmate::checkCharacter(trt)
  if (!isTRUE(good_trt)) {
    return(good_trt)
  }

  tau <- determine_tau(outcome, trt)
  if (length(trt) != 1 && length(trt) != tau) {
    return(paste0("'trt' should be of length 1 or ", tau))
  }

  good_cens <- checkmate::checkCharacter(cens, len = tau, null.ok = !any(is.na(x[[outcome]])))
  if (!isTRUE(good_cens)) {
    return(good_cens)
  }

  good_baseline <- checkmate::checkCharacter(baseline, null.ok = TRUE)
  if (!isTRUE(good_baseline)) {
    return(good_baseline)
  }

  good_time_vary <- checkmate::checkList(time_vary, types = c("NULL", "character"), len = tau, null.ok = TRUE)
  if (!isTRUE(good_time_vary)) {
    return(good_time_vary)
  }

  good_id <- checkmate::checkCharacter(id, len = 1, null.ok = TRUE)
  if (!isTRUE(good_id)) {
    return(good_id)
  }

  vars_exist <- checkmate::checkSubset(c(trt, outcome, baseline, unlist(time_vary), cens, id), names(x))
  if (!isTRUE(vars_exist)) {
    return(vars_exist)
  }

  for (t in 1:tau) {
    ci <- censored(x, cens, t)$j
    di <- at_risk(x, cens, t)
    trt_t <- ifelse(length(trt) == 1, trt, trt[t])
    data_t <- x[ci & !di, c(trt_t, baseline, unlist(time_vary[t])), drop = FALSE]

    if (any(is.na(data_t))) {
      return("Missing data found in treatment and/or covariate nodes for uncensored observations.")
    }
  }

  TRUE
}

assert_lmtp_data <- assertLmtpData <- checkmate::makeAssertionFunction(check_lmtp_data)

check_reserved_names <- function(x) {
  bad_names <- c("lmtp_id", "tmp_lmtp_stack_indicator", "tmp_lmtp_scaled_outcome") %in% names(x)
  if (!any(bad_names)) {
    return(TRUE)
  }
  "'lmtp_id', 'tmp_lmtp_stack_indicator', and 'tmp_lmtp_scaled_outcome' are reserved variable names."
}

assert_reserved_names <- assertReservedNames <- checkmate::makeAssertionFunction(check_reserved_names)

fix_censoring_ind <- function(data, C = NULL, tau) {
  if (!is.null(C)) {
    out <- as.list(data)
    for (t in 1:tau) {
      out[[C[[t]]]] <- ifelse(is.na(out[[C[[t]]]]), 0, out[[C[[t]]]])
    }
    return(as.data.frame(out))
  }

  data
}

check_variation <- function(outcome, learners) {
  if (sd(outcome) < .Machine$double.eps) {
    return(sl3::Lrnr_mean$new())
  }
  learners
}

check_outcome_type <- function(fits, ref, type) {
  if (type == "additive") {
    check <- TRUE
  }

  if (type %in% c("rr", "or")) {

    if (is.lmtp(ref)) {
      fits[["ref"]] <- ref
    }

    check <- all(lapply(fits, function(x) x[["outcome_type"]]) == "binomial")
  }

  if (isFALSE(check)) {
    stop(toupper(type), " contrast specified but one or more outcome types are non-binary.",
         call. = FALSE)
  }
}

check_lmtp_type <- function(fits, ref) {
  if (is.lmtp(ref)) {
    fits[["ref"]] <- ref
  }

  types <- lapply(fits, function(x) x[["estimator"]])
  check <- all(types %in% c("TMLE", "SDR"))

  if (isFALSE(check)) {
    stop("Contrasts not implemented for substitution/IPW estimators.",
         call. = FALSE)
  }
}

check_ref_type <- function(ref, type) {
  if (!is.lmtp(ref)) {
    if (class(ref) %in% c("numeric", "integer", "double")) {
      if (length(ref) != 1) {
        stop("Reference value should be a single object.",
             call. = FALSE)
      }

      message("Non-estimated reference value, defaulting type = 'additive'.")
      return("additive")
    }

    stop("Reference must either be a single numeric value or another lmtp object.",
         call. = FALSE)
  }

  out <- type
  out
}

check_trt_length <- function(trt, time_vary = NULL, cens = NULL, tau) {
  if (length(trt) == tau) {
    if (!is.null(time_vary) & !is.null(cens)) {
      if (tau == length(time_vary) & tau == length(cens)) {
        set_lmtp_options("trt", "standard")
        return(trt)
      }

      stop("It appears there is a mismatch in your data. Make sure parameters `trt`, `time_vary`, and `cens` are of the same length.",
           call. = FALSE)
    }

    if (!is.null(time_vary) & is.null(cens)) {
      if (tau == length(time_vary)) {
        set_lmtp_options("trt", "standard")
        return(trt)
      }

      stop("It appears there is a mismatch in your data. Make sure parameters `trt` and `time_vary` are of the same length.",
           call. = FALSE)
    }

    if (is.null(time_vary & !is.null(cens))) {
      if (tau == length(cens)) {
        set_lmtp_options("trt", "standard")
        return(trt)
      }

      stop("It appears there is a mismatch in your data. Make sure parameters `trt` and `cens` are of the same length",
           call. = FALSE)
    }

    set_lmtp_options("trt", "standard")
    return(trt)
  }

  set_lmtp_options("trt", "point.wise")
  rep(trt, tau)
}

check_at_risk <- function(outcomes, tau) {
  if (length(outcomes) == 1) {
    return(NULL)
  }

  if (length(outcomes) == tau) {
    return(outcomes[1:tau - 1])
  }

  stop("It appears there is a mismatch between the length of `outcome` and other parameters.",
       call. = FALSE)
}

 check_glm_outcome <- function(outcome_type) {
   if (is.null(outcome_type)) {
     return("gaussian")
   }

   if (outcome_type == "continuous") {
     return("gaussian")
   }

   if (outcome_type == "binomial") {
     return(outcome_type)
   }
 }

check_mult_outcomes <- function(outcome, outcome_type) {
  if (outcome_type == "survival") {
    # if (length(outcome) == 1) {
    #   stop("'outcome_type' set to survival, but the length of 'outcome' is one.", call. = FALSE)
    # }
    return()
  }

  if (length(outcome) > 1) {
    stop("The length of 'outcome' is greater than one, but 'outcome_type' not set to survival.", call. = FALSE)
  }
}

check_is_binary <- function(data, outcome, outcome_type) {
  if (outcome_type %in% c("binomial", "survival")) {
    vals <- lapply(outcome, function(x) {
      as.character(unique(na.omit(data[[x]])))
    })

    if (!all(unlist(lapply(vals, function(x) all(x %in% c("0", "1")))))) {
      stop("Only 0 and 1 alllowed in outcome variables if 'outcome_type' set to binomial or survival.",
           call. = FALSE)
    }
  }
}

check_factors <- function(data, trt, baseline, nodes) {
  vars <- c(trt, baseline, unlist(nodes))
  if (any(unlist(lapply(data[, vars], is.factor)))) {
    warning("Some of your variables appear to be factors. Make sure your SuperLearner library is capable of handling factors!",
            call. = FALSE)
  }
}

check_shifted <- function(data, shifted, outcome, baseline, nodes, cens, survival = FALSE) {
  unchngd <- c(outcome, baseline, unlist(nodes))
  shifted <- as.data.frame(shifted)

  if (survival) {
    for (outcomes in outcome) {
      data.table::set(shifted, j = outcomes, value = convert_to_surv(shifted[[outcomes]]))
    }
  }

  if (!(identical(data[unchngd], shifted[unchngd]))) {
    stop("If supplying data to `shifted`, the only columns that can be different between `data` and `shifted` are those indicated in `trt` and `cens`", .call = FALSE)
  }

  if (is.null(cens)) {
    return(shifted)
  }

  for (i in cens) {
    if (!(all(shifted[[i]] == 1))) {
      stop("If supplying data to `shifted`, censoring columns should be set to 1.", call. = FALSE)
    }
  }

  shifted
}

