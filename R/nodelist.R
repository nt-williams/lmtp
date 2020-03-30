
create_node_list <- function(A, history) {

  out <- list()
  tau <- length(A)
  A <- Reduce(c, A, accumulate = T)
  for (i in 1:tau) {
      out[[i]] <- c(history[[i]], A[[i]])
  }

  # returns
  out
}



