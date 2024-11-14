pivot <- function(x, Vars, tau) {
  vars <- lapply(1:tau, \(t) c("..i..lmtp_id", Vars$time(t)))
  to_bind <- lapply(vars, \(vars) x[, vars])
  to_bind <- lapply(to_bind, \(data) setNames(data, Vars$rename(names(data))))
  to_bind <- lapply(1:tau, function(t) {
    to_bind[[t]]$..i..wide_id <- 1:nrow(to_bind[[t]])
    to_bind[[t]]$time <- as.character(t)
    to_bind[[t]]
  })
  # NEED TO THINK ABOUT HOW TO HANDLE THE EDGE OF CASE OF NON-SURVIVAL TIME-VARYING Y
  longer <- do.call(rbind, to_bind)
  # longer <- longer[do.call(what = order, args = longer[,c("..i..wide_id", "time")]), ]

  # Create a new column with lagged Y
  if (!is.null(Vars$N)) {
    longer$..i..N <- ave(longer$..i..Y_1, longer$..i..wide_id, FUN = \(x) c(1, x[-length(x)]))
  } else {
    longer$..i..N <- rep(1, nrow(longer))
  }

  row.names(longer) <- NULL
  longer
}
