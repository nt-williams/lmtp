
#' Transform Survival Data Into the Correct Format
#'
#' @param data
#' @param trt
#' @param status
#' @param time
#' @param horizon
#'
#' @return
#' @export
#'
#' @importFrom data.table data.table `:=` fifelse .I dcast
#'
#' @examples
prep_survival_data <- function(formula, data, target, horizon = NULL) {
  time <- get_time(formula)
  status <- get_status(formula)
  covar <- get_covar(formula, target)
  nobs <- nrow(data)
  max_time <- horizon
  id <- 1:nobs

  if (is.null(max_time)) {
    max_time <- max(data[[time]])
  }

  all_time <- rep(1:max_time, nobs)
  evnt <- cens <- rep(NA, nobs*max_time)
  risk_evnt <- risk_cens <- 1*(all_time == 1)

  for (t in 1:max_time) {
    cens[all_time == t] <- (1 - data[[status]]) * (data[[time]] == t)
    evnt[all_time == t] <- data[[status]] * (data[[time]] == t)
  }

  long <-
    data.table(id = rep(id, each = max_time),
               data[as.numeric(gl(nobs, max_time)), ],
               all_time = all_time,
               outcome = evnt, status = cens)

  fevnt <- long[, .I[(.I[which(outcome == 1)]) < .I], by = id]$V1
  fcens <- long[, .I[(.I[which(status == 1)]) < .I], by = id]$V1

  long[fevnt, outcome := 1]
  long[fcens, `:=`(outcome = NA, status = 1)]
  long[, status := fifelse(status == 1, 0, 1)]
  long[all_time == max(all_time) & !is.na(outcome), status := 1]

  form <- paste(paste(c("id", trt, covar), collapse = "+"), "~ all_time")
  wide <- dcast(long, form, value.var = c("status", "outcome"))
  dlte <- c("outcome_1", paste0("status_", max_time))

  list(data = wide[, !dlte, with = FALSE],
       trt = target,
       baseline = covar,
       outcome = paste0("outcome_", 2:max_time),
       cens = paste0("status_", 1:(max_time - 1)))
}

get_status <- function(formula) {
  all.vars(formula[[2]])[[2]]
}

get_time <- function(formula) {
  all.vars(formula[[2]])[[1]]
}

get_covar <- function(formula, target) {
  setdiff(all.vars(formula[[3]]), target)
}
