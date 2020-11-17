
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
prep_survival_data <- function(formula, data, target, horizon = NULL, id = NULL) {
  time1 <- get_time1(formula)
  time2 <- get_time2(formula)
  status <- get_status(formula)
  covar <- get_covar(formula, target)
  max_time <- horizon

  setDT(data)

  if (is.null(id)) {
    data[, id := 1:nrow(data)]
    nobs <- nrow(data)
  } else {
    nobs <- length(unique(data[[id]]))
  }

  all_time <- rep(1:max_time, nobs)
  evnt <- cens <- rep(NA, nobs*max_time)

  browser()

  # tidyr::complete(dplyr::group_by(data, id), start = tidyr::full_seq(0:max(stop), 1)) %>%
  #   tidyr::fill(everything(), .direction = "downup")

  filled <-
    setnames(do.call(CJ, list(data[[id]], 1:max_time, unique = TRUE)), c(id, time2))
  long <- data[filled, on = c(id, time2)]

  long[, (target) := nafill(get(target), "nocb"), "id"]

  set(long, j = target, value = nafill(long[[target]], "nocb"))

  # complete_time(data, list(id = id, stop = 0:max_time))
  #
  # data[CJ(id = id, stop = 0:max(stop), unique = TRUE),
  #      on = .(id, stop)
  #      ][, trt := nafill(trt, "nocb")]

  lt <- data[, .I[.N], by = id]$V1

  for (t in 1:max_time) {
    cens[all_time == t] <- (1 - data[lt, status]) * (data[lt, get(time2)] == t)
    evnt[all_time == t] <- data[lt, status] * (data[lt, get(time2)] == t)
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
  vars <- all.vars(formula[[2]])
  if (length(vars) > 2) {
    return(vars[[3]])
  }
  vars[[2]]
}

get_time1 <- function(formula) {
  all.vars(formula[[2]])[[1]]
}

get_time2 <- function(formula) {
  vars <- all.vars(formula[[2]])
  if (length(vars) == 2) {
    return(NULL)
  }
  vars[[2]]
}

get_covar <- function(formula, target) {
  setdiff(all.vars(formula[[3]]), target)
}
