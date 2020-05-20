
check_for_sl3 <- function(test = FALSE) {
  tryCatch(
    if (isTRUE(test)) {
      stop()
    } else {
      has_sl3 <- "sl3" %in% rownames(utils::installed.packages())
      if (isFALSE(has_sl3)) stop()
      else on.exit()
    }, error = function(e) {
      no_sl3()
    }
  )
}

no_stderr_warning <- function(estimator) {
  cat("\n")
  cli::cli_alert_warning("Standard errors aren't provided for the {estimator} estimator.")
}

no_sl3 <- function() {
  cli::cli_text("Remote package, {.pkg sl3}, not detected.")
  cli::cli_text("{.pkg sl3} can be installed with: {.code remotes::install_github('tlverse/sl3')}")
}

check_sd <- function(x, learner_stack) {
  if (sd(x) > .Machine$double.eps) {
    out <- learner_stack
  } else {
    out <- sl3::Lrnr_mean
  }

  return(out)
}

check_censoring <- function(data, C, Y) {
  if (any(is.na(data[[Y]])) & is.null(C)) {
    stop("Missing outcomes detected and censoring nodes not indicated.", call. = FALSE)
  }
}

check_missing_data <- function(data, trt, nodes, baseline, cens, tau) {
  for (t in 1:tau) {
    i <- create_censoring_indicators(data, cens, t)$i
    if (any(is.na(as.matrix(data[i, c(check_trt_length(trt, tau)[t], baseline, unlist(nodes[t]))])))) {
      browser()
      stop("Missing data found in treatment and/or covariate nodes. Either impute (recommended) or only use observations with complete treatment and covariate data.",
           call. = F)
    }
  }
}

check_for_variables <- function(data, trt, outcome, baseline, nodes, cens) {
  vars <- c(trt, outcome, baseline, unlist(nodes), cens)
  if (!all(vars %in% names(data))) {
    warn <- vars[which(!(vars %in% names(data)))]
    stop("Variable(s) ", paste(warn, collapse = ", "), " not found in data.",
         call. = F)
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

  return(data)
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

check_extreme_ratio <- function(ratio) {
  return(apply(ratio, 2, function(x) pmin(x, quantile(x, 0.999))))
}

check_variation <- function(data, outcome, learners) {
  if (sd(data[, outcome]) < .Machine$double.eps) {
    learners <- sl3::make_learner(sl3::Lrnr_mean)
  }
  return(learners)
}

check_outcome_type <- function(fits, ref, type) {
  if (type == "additive") {
    check <- TRUE
  } else if (type %in% c("rr", "or")) {

    if (is.lmtp(ref)) {
      fits[["ref"]] <- ref
    }

    types <- lapply(fits, function(x) x[["outcome_type"]])
    check <- all(types == "binomial")
  }

  if (isFALSE(check)) {
    stop(toupper(type), " contrast specified but one or more outcome types are non-binary.",
         call. = F)
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
         call. = F)
  }

}

check_ref_type <- function(ref, type) {
  if (!is.lmtp(ref)) {
    if (class(ref) %in% c("numeric", "integer", "double")) {
      if (length(ref) != 1) {
        stop("Reference value should be a single object.",
             call. = F)
      }
      message("Non-estimated reference value, defaulting type = 'additive'.")
      out <- "additive"
    } else {
      stop("Reference must either be a single numeric value or another lmtp object.",
           call. = F)
    }
  } else {
    out <- type
  }

  return(out)
}

check_trt_length <- function(trt, tau) {
  if (length(trt) == tau) {
    set_lmtp_options("trt", "standard")
    return(trt)
  } else if (length(trt) == 1) {
    set_lmtp_options("trt", "point.wise")
    return(rep(trt, tau))
  } else {
    stop("Treatment nodes should either be the same length as nodes, or of length 1.",
         call. = F)
  }
}

check_deterministic <- function(outcomes, tau) {
  if (length(outcomes) == 1) {
    return(NULL)
  } else if (length(outcomes) == tau + 1) {
    return(outcomes[1:tau])
  } else {
    stop("Outcome argument must be of length 1, or in the case of survival analyis, of length tau + 1, with nodes 1 through tau set to the intermediate outcomes.",
         call. = F)
  }
}

check_folds <- function(V) {
  if (V > 1) {
    on.exit()
  } else {
    stop("The number of folds must be greater than 1.", call. = F)
  }
}
