#' Time To Event Last Outcome Carried Forward
#'
#' A helper function to prepare survival data for use with LMTP estimators
#' by imputing outcome nodes using last outcome carried forward when an observation
#' experiences the event before the end-of-follow-up.
#'
#' @param data The dataset to modify.
#' @param outcomes A vector of outcome nodes ordered by time.
#'
#' @return A modified dataset with future outcome nodes set to 1 if an observation
#'  experienced an event at any previous time point.
#'
#' @importFrom data.table as.data.table `:=` .SD
#' @export
#' @examples
#' event_locf(sim_point_surv, paste0("Y.", 1:6))
event_locf <- function(data, outcomes) {
  DT <- as.data.table(data)
  tau <- length(outcomes)
  for (j in outcomes[1:(tau - 1)]) {
    modify <- setdiff(outcomes[match(j, outcomes):tau], j)
    DT[get(j) == 1 & !is.na(get(j)), (modify) := lapply(.SD, function(x) 1), .SDcols = modify]
  }
  DT[]
  DT
}
