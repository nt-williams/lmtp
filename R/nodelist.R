
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

