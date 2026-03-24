# one_time_delay <- function(data, sequence) {
#   time <- length(sequence)
#   if (time == 1) return(rep(0, nrow(data)))
#
#   current <- data[[sequence[time]]]
#   previous <- data[[sequence[time - 1]]]
#
#   current * previous + (1 - current) * current
# }

one_time_delay <- function(data, sequence) {
  time <- length(sequence)
  if (time == 1) return(rep(0, nrow(data)))

  at <- data[[sequence[time]]]
  abartm1 <- data[, sequence[seq_len(time - 1)], drop = FALSE]
  test <- at * as.numeric(rowSums(abartm1) == 0)
  ifelse(test, 0, at)
}
