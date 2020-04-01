
create_node_list <- function(A, nodes, k = NULL) {

  out <- list()
  tau <- length(A)
  if (is.null(k) || k == 0) k <- Inf

  for (i in 1:tau) {
      out[[i]] <- c(nodes[[i]], A[i])
  }

  out <- paste(lapply(out, function(x) paste(x, collapse = ",")))
  out <- slider::slide(out, ~ .x, .before = k)
  out <- sapply(out, function(x) unlist(strsplit(x, ",")))

  # returns
  out
}

