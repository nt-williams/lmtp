LmtpWideTask <- R6Class("LmtpWideTask",
  inherit = LmtpTask,
  public = list(
    initialize = function(data, shifted, A, Y, W, L, C, id, weights, outcome_type, mtp, folds = 1) {
      super$initialize("wide", outcome_type, mtp)
      self$tau <- determine_tau(Y, A)

      self$col_roles <- LmtpWideVars$new(
        W = W,
        L = L,
        A = A,
        C = C,
        Y = Y,
        id = id,
        weights = weights,
        survival = outcome_type == "survival",
        tau = self$tau
      )

      self$backend <- private$as_lmtp_data(data)
      self$shifted <- private$as_lmtp_data(shifted)

      self$col_roles$id <- "lmtp_id"

      self$folds <- private$make_folds(folds)

      private$.row_copy <- 1:nrow(self$backend)
      private$.col_copy <- names(self$backend)
      self$active_rows <- private$.row_copy
      self$active_cols <- private$.col_copy
    },

    obs = function(t) {
      if (is.null(self$col_roles$C) | t == 0) {
        return(invisible(self))
      }

      var <- self$col_roles$C[t]

      obs <- self$backend[[var]][self$active_rows]
      self$active_rows(intersect(self$active_rows, which(obs == 1)))
      invisible(self)
    },

    at_risk = function(t) {
      if (!(self$outcome_type == "survival") || t == 1) {
        return(invisible(self))
      }

      risk <- self$col_roles$Y[1:(length(self$col_roles$Y) - 1)]
      risk <- self$backend[[risk[t - 1]]][private$.row_roles, ]

      self$active_rows(intersect(private$.row_roles, which(risk == 1 & !is.na(risk))))

      invisible(self)
    },

    followed_rule = function(t) {
      if (self$mtp) {
        return(rep(TRUE, self$nrow()))
      }

      if (length(self$col_roles$A) > 1) {
        trt_t <- self$col_roles$A[[t]]
      } else {
        trt_t <- self$col_roles$A[[1]]
      }

      x <- self$select(trt_t)$data(reset = FALSE)
      y <- self$select(trt_t)$data("shifted", reset = FALSE)

      mapply(function(x, y) isTRUE(all.equal(x, y)), as.list(x), as.list(y))
    },

    shift = function(t) {
      shifted <- self$data(reset = FALSE)
      A <- self$col_roles$A

      if (length(A) > 1) {
        At <- A[[t]]
      } else {
        At <- A[[1]]
      }

      shifted[, At] <- self$data("shifted", reset = FALSE)[, At]
      shifted
    },

    stack = function(t) {
      shifted_half <- natural <- self$data(reset = FALSE)

      A <- self$col_roles$A
      C <- self$col_roles$C

      if (length(A) > 1 || t == 1) {
        shifted_half[, A[[t]]] <- self$data("shifted", reset = FALSE)[, A[[t]]]
      }

      if (!is.null(C)) {
        shifted_half[[C[t]]] <- self$data("shifted", reset = FALSE)[[C[t]]]
      }

      out <- rbind(natural, shifted_half)
      out[["tmp_lmtp_stack_indicator"]] <- rep(c(0, 1), each = nrow(natural))
      out
    }
  ),
  private = list(
    as_lmtp_data = function(data) {
      assert_subset(self$col_roles$all(), names(data))

      assert_lmtp_data(
        data,
        self$col_roles$A,
        self$col_roles$Y,
        self$col_roles$W,
        self$col_roles$L,
        self$col_roles$C,
        self$col_roles$id
      )

      data <- data.table::copy(as.data.frame(data))
      data$lmtp_id <- create_ids(data, self$col_roles$id)
      data <- fix_censoring_ind(data, self$col_roles$C)

      if (self$outcome_type == "survival") {
        for (y in self$col_roles$Y) {
          data.table::set(data, j = y, value = convert_to_surv(data[[y]]))
        }
      }

      Y_tau <- self$col_roles$Y[length(self$col_roles$Y)]
      private$bounds <- y_bounds(data[[Y_tau]], self$outcome_type)
      data$tmp_lmtp_scaled_outcome <- scale_y(data[[Y_tau]], private$bounds)

      if (!is.null(self$col_roles$weights)) {
        wts <- data[[self$col_roles$weights]]
        if (!is_normalized(wts)) {
          data[[self$col_roles$weights]] <- wts / mean(wts)
        }
      }

      data
    }
  )
)

LmtpWideTaskSplit <- R6Class("LmtpWideTaskSplit",
  inherit = LmtpWideTask,
  public = list(
    type = NULL,
    initialize = function(x, type) {
      self$type <- type
      self$backend <- x$data(reset = FALSE)
      self$shifted <- x$data("shifted")
      self$outcome_type <- x$outcome_type
      self$col_roles <- x$col_roles
      self$tau <- x$tau
      self$mtp <- x$mtp
      private$.row_copy <- 1:nrow(self$backend)
      private$.col_copy <- names(self$backend)
      self$active_rows <- private$.row_copy
      self$active_cols <- private$.col_copy
    }
  )
)
