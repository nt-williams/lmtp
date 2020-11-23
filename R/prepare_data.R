# #' Transform Survival Data Into the Correct Format
# #'
# #' @param data
# #' @param trt
# #' @param status
# #' @param time
# #' @param horizon
# #'
# #' @return
# #'
# #' @importFrom data.table data.table `:=` fifelse .I dcast
# #'
# #' @examples
# prep_survival_data <- function(formula, data, target, horizon = NULL, id = NULL) {
#   setDT(data)
#   time1 <- get_time1(formula)
#   time2 <- get_time2(formula)
#   status <- get_status(formula)
#   covar <- get_covar(formula, target)
#
#   if (!is.null(time2)) {
#     time <- time2
#     min_time <- min(data[[time1]])
#   } else {
#     time <- time1
#     min_time <- 1
#   }
#
#   max_time <- max(data[[time]])
#   time_seq <- min_time:max_time
#   all_vars <- c(id, time1, time2, target, status, covar)
#   data <- data[, .SD, .SDcols = all_vars]
#
#   if (is.null(id)) {
#     data[, id := 1:nrow(data)]
#     nobs <- nrow(data)
#     id <- "id"
#   } else {
#     nobs <- length(unique(data[[id]]))
#   }
#
#   all_time <- rep(time_seq, nobs)
#   evnt <- cens <- rep(NA, nobs*max_time)
#
#   filled <-
#     setnames(do.call(CJ, list(data[[id]], time_seq, unique = TRUE)), c(id, time))
#   long <- data[filled, on = c(id, time)]
#
#   if (is.null(time2)) {
#     fx <- c(target, covar)
#     long[, (fx) := lapply(.SD, zoo::na.locf, na.rm = FALSE, fromLast = TRUE), by = id, .SDcols = fx
#          ][, (fx) := lapply(.SD, zoo::na.locf, na.rm = FALSE), by = id, .SDcols = fx]
#   }
#
#   if (!is.null(time2)) {
#     static <- unlist(lapply(covar, is_baseline, data = data, id = id))
#     cnst <- covar[static]
#     time_vary <- covar[!static]
#     vars <- c(target, covar)
#     long[, (vars) := lapply(.SD, zoo::na.locf, na.rm = FALSE, fromLast = TRUE), by = id, .SDcols = vars
#          ][, (cnst) := lapply(.SD, zoo::na.locf, na.rm = FALSE), by = id, .SDcols = cnst]
#   }
#
#   lt <- data[, .I[.N], by = id]$V1
#   for (t in time_seq) {
#     cens[all_time == t] <- (1 - data[lt][[status]]) * (data[lt][[time]] == t)
#     evnt[all_time == t] <- data[lt][[status]] * (data[lt][[time]] == t)
#   }
#
#   long[, `:=`(outcome = evnt,
#               status = cens,
#               all_time = long[[time]])]
#
#   fevnt <- long[, .I[(.I[which(outcome == 1)]) <= .I], by = id]$V1
#   fcens <- long[, .I[(.I[which(status == 1)]) < .I], by = id]$V1
#
#   long[fevnt, `:=`(outcome = 1, status = NA)]
#   long[fcens, `:=`(outcome = NA, status = 1)]
#   long[, status := fifelse(status == 1, 0, 1)]
#
#   if (is.null(time2)) {
#     form <- paste(paste(c("id", target, covar), collapse = "+"), "~ all_time")
#     wide <- dcast(long[all_time <= horizon, ], form, value.var = c("status", "outcome"))
#     dlte <- c("outcome_1", paste0("status_", horizon))
#
#     out <- list(data = wide[, !dlte, with = FALSE],
#                 trt = target,
#                 baseline = covar,
#                 outcome = paste0("outcome_", (min_time + 1):horizon),
#                 cens = paste0("status_", min_time:(horizon - 1)))
#   }
#
#   if (!is.null(time2)) {
#     # long[all_time == max(all_time) & !is.na(outcome), status := 1]
#     form <- paste(paste(c(id, cnst), collapse = "+"), "~ all_time")
#     wide <- dcast(long[all_time <= horizon, ], form, value.var = c(target, time_vary, "status", "outcome"))
#     dlte <- c("outcome_1", paste0("status_", horizon))
#
#     out <- list(data = wide[, !dlte, with = FALSE],
#                 trt = paste0(target, "_", min_time:(horizon - 1)),
#                 baseline = cnst,
#                 time_vary = NULL,
#                 outcome = paste0("outcome_", 2:horizon),
#                 cens = paste0("status_", 1:(horizon - 1)))
#   }
#   return(out)
# }
#
# get_status <- function(formula) {
#   vars <- all.vars(formula[[2]])
#   if (length(vars) > 2) {
#     return(vars[[3]])
#   }
#   vars[[2]]
# }
#
# get_time1 <- function(formula) {
#   all.vars(formula[[2]])[[1]]
# }
#
# get_time2 <- function(formula) {
#   vars <- all.vars(formula[[2]])
#   if (length(vars) == 2) {
#     return(NULL)
#   }
#   vars[[2]]
# }
#
# get_covar <- function(formula, target) {
#   setdiff(all.vars(formula[[3]]), target)
# }
#
# is_baseline <- function(data, x, id) {
#   chck <- data[, length(unique(get(x))), by = id]$V1
#   if (all(chck == 1)) {
#     return(TRUE)
#   }
#   return(FALSE)
# }
#
# # minus_one <- function(var, x) {
# #   paste0(var, "_", as.numeric(sub("\\D+", "", x)) - 1)
# # }
