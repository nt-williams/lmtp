LmtpWideVars <- R6Class("LmtpWideVars",
  public = list(
    W = NULL,
    L = NULL,
    A = NULL,
    C = NULL,
    Y = NULL,
    id = NULL,
    weights = NULL,
    initialize = function(W, L, A, C, Y, id, weights, survival, tau) {
      assert_trt(A, tau)
      assert_character(W, null.ok = TRUE)
      assert_character(C, len = tau, null.ok = TRUE)
      assert_character(Y, len = if (survival) 1, min.len = if (survival) 2)
      assert_list(L, types = c("NULL", "character"), len = tau, null.ok = TRUE)
      assert_character(id, len = 1, null.ok = TRUE)
      assert_character(weights, len = 1, null.ok = TRUE)

      self$W <- W
      self$L <- L
      self$A <- A
      self$C <- C
      self$Y <- Y
      self$id <- id
      self$weights <- weights
    },

    #' Get all parent nodes for a variable
    history = function(var = c("L", "A", "Y"), t) {
      switch(
        match.arg(var),
        L = private$parents_L(t),
        A = private$parents_A(t),
        Y = private$parents_Y()
      )
    },

    all = function() {
      c(self$W, unlist(self$L), unlist(self$A), self$C, self$Y, self$id, self$weights, self$tmp)
    }
  ),
  private = list(
    parents_L = function(t) {
      if (t == 1) {
        return(self$W)
      }
      c(private$parents_A(t - 1), unlist(self$A[t - 1]))
    },
    parents_A = function(t) {
      c(private$parents_L(t), unlist(self$L[[t]]))
    },
    parents_Y = function() {
      c(self$W, unlist(self$L), unlist(self$A))
    }
  )
)
