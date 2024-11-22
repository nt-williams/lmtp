pivot <- function(x, Vars) {
  vars <- lapply(1:Vars$tau, \(t) c("..i..lmtp_id", Vars$time(t)))
  to_bind <- lapply(vars, \(vars) x[, vars])
  to_bind <- lapply(to_bind, \(data) setNames(data, Vars$rename(names(data))))
  to_bind <- lapply(1:Vars$tau, function(t) {
    to_bind[[t]]$..i..wide_id <- 1:nrow(to_bind[[t]])
    to_bind[[t]]$time <- factor(t)
    to_bind[[t]]
  })
  # NEED TO THINK ABOUT HOW TO HANDLE THE EDGE OF CASE OF NON-SURVIVAL TIME-VARYING Y
  longer <- do.call(rbind, to_bind)

  # Create a new column with lagged Y
  if (!is.null(Vars$N)) {
    longer$..i..N <- ave(longer$..i..Y_1, longer$..i..wide_id, FUN = \(x) c(1, x[-length(x)]))
  } else {
    longer$..i..N <- rep(1, nrow(longer))
  }

  if (Vars$tau > 1) {
    k <- min(Vars$k, Vars$tau)
    longer <- as.data.table(longer)

    to_lag <- grep("^(..i..L)|(..i..A)", names(longer), value = TRUE)
    for (l in 1:(k-1)) {
      if (k > 0) {
        newcols <- paste0(to_lag, "_lag", l)
        longer[, (newcols) := .lag(.SD, l), by = ..i..wide_id, .SDcols = to_lag]
      }
    }
  }

  if (is.null(Vars$C)) {
    longer$..i..C_1 <- rep(1, nrow(longer))
    longer$..i..C_1_lag <- rep(1, nrow(longer))
  } else {
    longer$..i..C_1_lag <- ave(longer$..i..C_1, longer$..i..wide_id, FUN = \(x) c(1, x[-length(x)]))
  }

  row.names(longer) <- NULL
  as.data.frame(longer)
}

.lag <- function(x, n) {
  data.table::shift(x, n = n, type = "lag", fill = 1L)
}
