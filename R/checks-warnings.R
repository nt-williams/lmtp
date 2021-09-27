no_stderr_warning <- function(estimator) {
  cat("\n")
  cli::cli_alert_warning("Standard errors aren't provided for the {estimator} estimator.")
}

check_censoring <- function(data, C, Y) {
  if (any(is.na(data[[Y]])) & is.null(C)) {
    stop("Missing outcomes detected and censoring nodes not indicated.", call. = FALSE)
  }
}

check_missing_data <- function(data, trt, outcome, time_vary, baseline, cens, tau) {
  for (t in 1:tau) {
    ci <- censored(data, cens, t)$j
    di <- at_risk(data, cens, t)
    if (any(is.na(as.matrix(data[ci & !di, c(check_trt_length(trt, time_vary, cens, tau)[t], baseline, unlist(time_vary[t]))])))) {
      stop("Missing data found in treatment and/or covariate nodes. Either impute (recommended) or only use observations with complete treatment and covariate data.",
           call. = FALSE)
    }
  }
}

check_for_variables <- function(data, trt, outcome, baseline, nodes, cens) {
  vars <- c(trt, outcome, baseline, unlist(nodes), cens)
  if (!all(vars %in% names(data))) {
    warn <- vars[which(!(vars %in% names(data)))]
    stop("Variable(s) ", paste(warn, collapse = ", "), " not found in data.",
         call. = FALSE)
  }
}

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

check_scaled_conflict <- function(data) {
  nn <- names(data)
  check <- "xyz" %in% nn

  tryCatch(
    if (check) {
      stop()
    } else {
      on.exit()
    },
    error = function(e) {
      stop(
        "A variable named `xyz` was detected in your data. This variable name is reserved for an interal `lmtp` process. Please rename this variable.",
        call. = FALSE
      )
    }
  )
}

check_variation <- function(outcome, learners) {
  if (sd(outcome) < .Machine$double.eps) {
    return("SL.mean")
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

check_folds <- function(V) {
  if (V <= 1) {
    stop("The number of folds must be greater than 1.", call. = FALSE)
  }
}

 check_time_vary <- function(time_vary = NULL) {
   if (!is.null(time_vary)) {
     if (!is.list(time_vary)) {
       stop("time_vary must be a list.", call. = FALSE)
     }
   }
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
    if (length(outcome) == 1) {
      stop("'outcome_type' set to survival, but the length of 'outcome' is one.", call. = FALSE)
    }
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

check_shifted <- function(data, shifted, outcome, baseline, nodes, cens) {
  unchngd <- c(outcome, baseline, unlist(nodes))
  shifted <- as.data.frame(shifted)

  if (!(all.equal(data[unchngd], shifted[unchngd]))) {
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

