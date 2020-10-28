
# for now just going to work with a simple case of
# a non time-varying treatment
prepare_data <- function(data, trt, status, time) {
  nobs <- nrow(data)
  max_time <- max(data[[time]])
  all_time <- rep(1:max_time, nobs)
  evnt <- cens <- rep(NA, nobs*max_time)
  risk_evnt <- risk_cens <- 1*(all_time == 1)

  for (t in 1:max_time) {
    cens[all_time == t] <- (1 - data[[status]]) * (data[[time]] == t)
    evnt[all_time == t] <- data[[status]] * (data[[time]] == t)
    risk_evnt[all_time == t] <- (data[[time]] >= t)
    risk_cens[all_time == t] <- (data[[time]] > t) * data[[status]] +
      (data[[time]] >= t) * (1 - data[[status]])
  }

  data.table::data.table(id = rep(1:nobs, each = max_time),
                         data[as.numeric(gl(nobs, max_time)), ],
                         all_time = all_time,
                         evnt, cens)
}
