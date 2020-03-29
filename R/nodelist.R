
create_node_list <- function(A, history) {

  out <- list()
  for (i in 1:length(history)) {
    for (j in 1:length(A)) {
      out[[i]] <- c(history[[i]], A[[j]])
    }
  }

  # returns
  out
}
