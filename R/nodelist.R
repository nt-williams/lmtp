
#' Create a node list specification
#'
#' Creates a node list specification that is used by the provided estimators.
#' \code{create_node_list()} is not explicitly called by the analyst, rather
#' it is provided so the analyst can confirm how estimators will use variables
#' before actually performing the estimation procedure.
#'
#' @param trt A vector of column names of treatment variables.
#' @param nodes A list of length tau with the column names for new nodes to
#'  be introduced at each time point. The list should be ordered following
#'  the time ordering of the model.
#' @param baseline An optional vector of columns names for baseline covariates to be
#'  included for adjustment at every timepoint.
#' @param k An integer specifying how many previous time points should be
#'  used for estimation at the given time point. Default is \code{Inf},
#'  all time points.
#'
#' @return A list the same length of the nodes parameter with the variables
#' to be used for estimation at that given time point.
#' @export
#' @examples
#' a <- c("A_1", "A_2", "A_3", "A_4")
#' bs <- c("W_1", "W_2")
#' nodes <- list(c("L_1"), c("L_2"), c("L_3"), c("L_4"))
#'
#' # assuming no Markov property
#' create_node_list(a, nodes, bs, k = Inf)
#'
#' # assuming a Markov property
#' create_node_list(a, nodes, bs, k = 1)
create_node_list <- function(trt, nodes, baseline = NULL, k = Inf) {

  out <- list()
  tau <- length(nodes)

  if (is.null(k)) {
    k <- Inf
  }

  if (tau == 1 & k == Inf) {
    k <- 0
  }

  if (length(trt) == tau) {
    for (i in 1:tau) {
      out[[i]] <- c(nodes[[i]], trt[i])
    }
  } else {
    for (i in 1:tau) {
      out[[i]] <- c(nodes[[i]], trt)
    }
  }

  out <- paste(lapply(out, function(x) paste(x, collapse = ",")))
  out <- slider::slide(out, ~ .x, .before = k)
  out <- sapply(out, function(x) {
    . <- strsplit(x, ",")
    if (k == 0) .
    else unique(unlist(.))
  })

  if (!is.null(baseline)) {
    for (i in 1:tau) {
      out[[i]] <- c(baseline, out[[i]])
    }
  }

  # returns
  out
}

