pivot <- function(x, Vars, tau) {
  vars <- lapply(1:tau, \(t) c("._lmtp_id", Vars$time(t)))
  to_bind <- lapply(vars, \(vars) x[, vars])
  to_bind <- lapply(to_bind, \(data) setNames(data, Vars$rename(names(data))))
  to_bind <- lapply(1:tau, function(t) {
    to_bind[[t]]$.wide_id <- 1:nrow(to_bind[[t]])
    to_bind[[t]]$time <- t
    to_bind[[t]]
  })
  longer <- do.call(rbind, to_bind)
  row.names(longer) <- NULL
  longer
}
