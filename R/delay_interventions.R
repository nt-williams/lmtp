one_time_delay <- function(data, sequence) {
  time <- length(sequence)
  if (time == 1) return(rep(0, nrow(data)))

  current <- data[[sequence[time]]]
  previous <- data[[sequence[time - 1]]]

  current * previous + (1 - current) * current
}

# # one_time_delay <- function(data, sequence) {
# #   time <- length(sequence)
# #   if (time == 1) return(rep(0, nrow(data)))
#
# #   at <- data[[sequence[time]]]
# #   if (time == 11) return(at)
# #   abartm1 <- data[, sequence[seq_len(time - 1)], drop = FALSE]
# #   test <- at * as.numeric(rowSums(abartm1) == 0)
# #   ifelse(test, 0, at)
# # }
#
# # # Really a two-time delay
# # one_time_delay <- function(data, sequence) {
# #   time <- length(sequence)
# #   if (time %in% c(1,2)) return(rep(0, nrow(data)))
#
# #   at <- data[[sequence[time]]]
# #   if (time == 11) return(at)
# #   abartm1 <- data[, sequence[seq_len(time - 2)], drop = FALSE]
# #   test <- at * as.numeric(rowSums(abartm1) == 0)
# #   ifelse(test, 0, at)
# # }
#
# # Really a delay to only the last time
# one_time_delay <- function(data, sequence) {
#   time <- length(sequence)
#   if (time %in% seq_len(10)) return(rep(0, nrow(data)))
#
#   at <- data[[sequence[time]]]
#   if (time == 11) return(at)
# }
