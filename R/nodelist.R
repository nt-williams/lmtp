
#' Create a node list specification
#'
#' @param A A vector of column names of treatment variables.
#' @param nodes A list of length tau with the column names for new nodes to
#'  be introduced at each time point. The list should be ordered following
#'  the time ordering of the model.
#' @param baseline An optional vector of columns names for baseline covariates to be
#'  included for adjustment at every timepoint. If \code{k = Inf}, should be \code{NULL}
#'  and these variables should be added to the first index of \code{nodes}.
#' @param k An integer specifying how many previous time points nodes should be
#'  used for estimation at the given time point. Default is \code{Inf},
#'  all time points.
#'
#' @return TODO
#' @export
create_node_list <- function(A, nodes, baseline = NULL, k = Inf) {

  out <- list()
  tau <- length(A)

  if (is.null(k)) {
    k <- Inf
  }

  if (tau == 1 & k == Inf) {
    k <- 0
  }

  if (!is.null(baseline) & k == Inf) {
    stop("Non-Markov process, set `baseline = NULL` and add baseline covariates to `nodes[[1]]`",
         call. = FALSE)
  }

  for (i in 1:tau) {
    out[[i]] <- c(nodes[[i]], A[i])
  }

  out <- paste(lapply(out, function(x) paste(x, collapse = ",")))
  out <- slider::slide(out, ~ .x, .before = k)
  out <- sapply(out, function(x) {
    . <- strsplit(x, ",")
    if (k == 0) .
    else unlist(.)
  })

  if (!is.null(baseline)) {
    for (i in 1:tau) {
      out[[i]] <- c(baseline, out[[i]])
    }
  }

  # returns
  out
}

