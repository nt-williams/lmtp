LmtpVars <- R6Class("LmtpVars",
  public = list(
    W = NULL,
    L = NULL,
    A = NULL,
    C = NULL,
    Y = NULL,
    N = NULL,
    initialize = function(W, L, A, C, Y, tau, k = Inf) {
      assert_trt(A, tau)
      assert_character(W, null.ok = TRUE)
      assert_character(C, len = tau, null.ok = TRUE)
      assert_character(Y, min.len = 1, max.len = tau)
      assert_list(L, types = c("NULL", "character"), len = tau, null.ok = TRUE)
      assert_number(k, lower = 0, upper = Inf)

      self$W <- W
      self$L <- L
      self$A <- A
      self$C <- C
      self$Y <- Y

      if (length(self$Y) > 1) {
        self$N <- self$Y[1:length(self$Y) - 1]
      }

      private$tau <- tau
      private$k <- k
    },

    history = function(var = c("L", "A"), t) {
      private$l <- t - private$k - 1
      private$.var <- match.arg(var)
      if (private$.var == "A" && private$tau > 1 && length(self$A) == 1) {
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
        if (x %in% unlist(self$A)) prefix <- "._A"
        else if (x %in% unlist(self$L)) prefix <- "._L"
        else if (x %in% self$C) prefix <- "._C"
        else if (x %in% self$Y) prefix <- "._Y"
        else if (x %in% self$W) return(self$W[which(self$W == x)])
        else return(x)

        if (prefix == "._L" | (prefix == "._A" && is.list(self$A))) {
          vars <- self[[gsub("\\._", "", prefix)]]
          suffix <- which(x == vars[[which(sapply(vars, \(vars) x %in% vars))]])
        } else {
          suffix <- 1
        }

        paste0(prefix, "_", suffix)
      })
    },

    time = function(t) {
      A <- unlist(self$A[t])
      Y <- self$Y[t]
      if (is.na(A)) A <- self$A[1]
      if (is.na(Y)) Y <- self$Y[1]
      c(self$W, unlist(self$L[t]), A, self$C[t], Y)
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
