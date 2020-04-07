
#' Create a node list specification
#'
#' @param A A vector of column names of treatment variables.
#' @param nodes A list of length tau with the column names for new nodes to
#'  be introduced at each time point. The list should be ordered following
#'  the time ordering of the model.
#' @param k An integer specifying how many previous time points nodes should be
#'  used for estimation at the given time point. Default is \code{Inf},
#'  all time points.
#'
#' @return TODO
#' @export
create_node_list <- function(A, nodes, k = Inf) {

  out <- list()
  tau <- length(A)

  if (is.null(k)) k <- Inf
  if (tau == 1 & k == Inf) k <- 0

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

  # returns
  out
}

