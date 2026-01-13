one_time_delay <- function(data, sequence) {
  time <- length(sequence)
  if (time == 1) return(rep(0, nrow(data)))

  current <- data[[sequence[time]]]
  previous <- data[[sequence[time - 1]]]

  current * previous + (1 - current) * current
}
