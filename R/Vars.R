LmtpWideVars <- R6Class("LmtpWideVars",
  public = list(
    W = NULL,
    L = NULL,
    A = NULL,
    C = NULL,
    Y = NULL,
    id = NULL,
    weights = NULL,
    initialize = function(W, L, A, C, Y, id, weights, survival, tau, k = Inf) {
      assert_trt(A, tau)
      assert_character(W, null.ok = TRUE)
      assert_character(C, len = tau, null.ok = TRUE)
      assert_character(Y, len = if (!survival) 1, min.len = if (survival) 2)
      assert_list(L, types = c("NULL", "character"), len = tau, null.ok = TRUE)
      assert_character(id, len = 1, null.ok = TRUE)
      assert_character(weights, len = 1, null.ok = TRUE)
      checkmate::assertNumber(k, lower = 0, upper = Inf)

      self$W <- W
      self$L <- L
      self$A <- A
      self$C <- C
      self$Y <- Y
      self$id <- id
      self$weights <- weights

      private$tau <- tau
      private$k <- k
      # private$l <- tau - (tau - k - 1)
    },

    history = function(var = c("L", "A"), t) {
      private$l <- t - private$k - 1
      private$.var <- match.arg(var)
      ans <- switch(private$.var,
        L = private$parents_L(t),
        A = private$parents_A(t)
      )
      private$.var <- NULL
      private$l <- NULL
      as.vector(na.omit(c(self$W, ans)))
    },

    all = function() {
      c(self$W, unlist(self$L), unlist(self$A), self$C, self$Y, self$id, self$weights)
    }
  ),
  private = list(
    tau = NULL,
    k = NULL,
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
