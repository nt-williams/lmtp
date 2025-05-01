LmtpVars <- R6Class("LmtpVars",
  public = list(
    W = NULL,
    L = NULL,
    A = NULL,
    C = NULL,
    D = NULL,
    Y = NULL,
    N = NULL,
    tau = NULL,
    k = NULL,
    initialize = function(W, L, A, C, D, Y, outcome_type, tau, k = Inf) {
      assert_trt(A, tau)
      assert_character(W, null.ok = TRUE)
      assert_character(C, len = tau, null.ok = TRUE)
      assert_character(D, len = tau, null.ok = TRUE) # competing risk indicators
      assert_character(Y, min.len = ifelse(outcome_type == "survival", 2, 1))
      assert_list(L, types = c("NULL", "character"), len = tau, null.ok = TRUE)
      assert_number(k, lower = 0, upper = Inf)

      self$W <- W
      self$L <- L
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
      private$l <- t - self$k - 1
      private$.var <- match.arg(var)
      if (private$.var == "A" && self$tau > 1 && length(self$A) == 1) {
        return(as.vector(na.omit(self$W)))
      }
      ans <- switch(private$.var,
        L = private$parents_L(t),
        A = private$parents_A(t)
      )
      private$.var <- NULL
      private$l <- NULL
      as.vector(na.omit(c(self$W, ans)))
    },

    all = function() {
      c(self$W, unlist(self$L), unlist(self$A), self$C, self$Y, self$N)
    },

    rename = function(x) {
      sapply(x, function(x) {
        if (x %in% unlist(self$A)) prefix <- "..i..A"
        else if (x %in% unlist(self$L)) prefix <- "..i..L"
        else if (x %in% self$C) prefix <- "..i..C"
        else if (x %in% self$D) prefix <- "..i..D"
        else if (x %in% c(self$N, self$Y)) prefix <- "..i..Y"
        else if (x %in% self$W) return(self$W[which(self$W == x)])
        else return(x)

        if (prefix == "..i..L" | (prefix == "..i..A" && is.list(self$A))) {
          vars <- self[[gsub("\\..i..", "", prefix)]]
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
      c(self$W, unlist(self$L[t]), A, self$C[t], self$D[t], Y)
    }
  ),
  private = list(
    l = NULL,
    .var = NULL,
    parents_L = function(t) {
      if (t == 1) {
        return(invisible())
      }

      if (t == private$l) {
        return(unlist(self$A[t - 1]))
      }

      c(private$parents_A(t - 1), unlist(self$A[t - 1]))
    },
    parents_A = function(t) {
      if (t == private$l) {
        if (private$.var == "L") {
          return(unlist(self$L[[t]]))
        }
        return(invisible())
      }

      c(private$parents_L(t), unlist(self$L[[t]]))
    }
  )
)
