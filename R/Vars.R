LmtpVars <- R6Class("LmtpVars",
  public = list(
    W = NULL,
    W_trt = NULL,
    W_cens = NULL,
    L = NULL,
    L_trt = NULL,
    L_cens = NULL,
    A = NULL,
    C = NULL,
    D = NULL,
    Y = NULL,
    N = NULL,
    tau = NULL,
    k = NULL,
    initialize = function(W, L, A, C, D, Y, outcome_type, tau, k = Inf) {
      assert_trt(A, tau)

      # Handle W (baseline) - accept old format (character) or new format (named list)
      if (is.list(W) && !is.null(names(W))) {
        assert_subset(names(W), c("trt", "cens", "outcome"))
        assert_character(W$trt, null.ok = TRUE)
        assert_character(W$cens, null.ok = TRUE)
        assert_character(W$outcome, null.ok = TRUE)
        self$W <- W$outcome
        self$W_trt <- W$trt
        self$W_cens <- W$cens
      } else {
        assert_character(W, null.ok = TRUE)
        self$W <- W
        self$W_trt <- W
        self$W_cens <- W
      }

      # Handle L (time_vary) - accept old format (list of length tau) or new format (named list of 3)
      if (is.list(L) && !is.null(L) && !is.null(names(L))) {
        assert_subset(names(L), c("trt", "cens", "outcome"))
        assert_list(L$trt, types = c("NULL", "character"), len = tau, null.ok = TRUE)
        assert_list(L$cens, types = c("NULL", "character"), len = tau, null.ok = TRUE)
        assert_list(L$outcome, types = c("NULL", "character"), len = tau, null.ok = TRUE)
        self$L <- L$outcome
        self$L_trt <- L$trt
        self$L_cens <- L$cens
      } else {
        assert_list(L, types = c("NULL", "character"), len = tau, null.ok = TRUE)
        self$L <- L
        self$L_trt <- L
        self$L_cens <- L
      }

      assert_character(C, len = tau, null.ok = TRUE)
      assert_character(D, len = tau, null.ok = TRUE) # competing risk indicators
      assert_character(Y, min.len = ifelse(outcome_type == "survival", 2, 1))
      assert_number(k, lower = 0, upper = Inf)

      self$A <- A
      self$C <- C

      if (outcome_type == "survival") {
        self$Y <- last(Y)
        self$N <- Y[1:length(Y) - 1]
        self$D <- D
      } else {
        self$Y <- Y
      }

      self$k <- k
      self$tau <- tau
    },

    history = function(var = c("L", "A"), t) {
      var <- match.arg(var)
      W_use <- if (var == "A") self$W_trt else self$W
      L_use <- if (var == "A") self$L_trt else self$L
      private$compute_history(var, t, W_use, L_use)
    },

    history_cens = function(t) {
      # Censoring model uses same history structure as outcome (var="L") but cens-specific W/L
      private$compute_history("L", t + 1, self$W_cens, self$L_cens)
    },

    all = function() {
      unique(c(self$W, self$W_trt, self$W_cens,
               unlist(self$L), unlist(self$L_trt), unlist(self$L_cens),
               unlist(self$A), self$C, self$Y, self$N))
    },

    rename = function(x) {
      all_L <- unique(c(unlist(self$L), unlist(self$L_trt), unlist(self$L_cens)))
      all_W <- unique(c(self$W, self$W_trt, self$W_cens))
      sapply(x, function(x) {
        if (x %in% unlist(self$A)) prefix <- "..i..A"
        else if (x %in% all_L) prefix <- "..i..L"
        else if (x %in% self$C) prefix <- "..i..C"
        else if (x %in% self$D) prefix <- "..i..D"
        else if (x %in% c(self$N, self$Y)) prefix <- "..i..Y"
        else if (x %in% all_W) return(all_W[which(all_W == x)])
        else return(x)

        # Use outcome L for position lookup (backward compat)
        L_for_rename <- if (!is.null(self$L)) self$L else if (!is.null(self$L_trt)) self$L_trt else self$L_cens
        if (prefix == "..i..L" | (prefix == "..i..A" && is.list(self$A))) {
          vars <- if (prefix == "..i..L") L_for_rename else self$A
          suffix <- which(x == vars[[which(sapply(vars, function(vars) x %in% vars))]])
        } else {
          suffix <- 1
        }

        paste0(prefix, "_", suffix)
      })
    },

    time = function(t) {
      A <- unlist(self$A[t])
      Y <- c(self$N, self$Y)[t]
      if (all(is.na(A))) A <- unlist(self$A[1])
      if (is.na(Y)) Y <- self$Y[1]
      unique(c(self$W, self$W_trt, self$W_cens,
               unlist(self$L[t]), unlist(self$L_trt[t]), unlist(self$L_cens[t]),
               A, self$C[t], self$D[t], Y))
    }
  ),
  private = list(
    compute_history = function(var, t, W, L) {
      private$l <- t - self$k - 1
      private$.var <- var
      if (var == "A" && self$tau > 1 && length(self$A) == 1) {
        private$l <- NULL
        private$.var <- NULL
        return(as.vector(na.omit(W)))
      }
      ans <- switch(var,
        L = private$parents_L(t, L),
        A = private$parents_A(t, L)
      )
      private$l <- NULL
      private$.var <- NULL
      as.vector(na.omit(c(W, ans)))
    },
    l = NULL,
    .var = NULL,
    parents_L = function(t, L) {
      if (t == 1) {
        return(invisible())
      }

      if (t == private$l) {
        return(unlist(self$A[t - 1]))
      }

      c(private$parents_A(t - 1, L), unlist(self$A[t - 1]))
    },
    parents_A = function(t, L) {
      if (t == private$l) {
        if (private$.var == "L") {
          return(unlist(L[[t]]))
        }
        return(invisible())
      }

      c(private$parents_L(t, L), unlist(L[[t]]))
    }
  )
)
