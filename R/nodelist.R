#' Create a node list specification
#'
#' Creates a node list specification that is used by the provided estimators.
#' \code{create_node_list()} is not explicitly called by the analyst, rather
#' it is provided so the analyst can confirm how estimators will use variables
#' before actually performing the estimation procedure.
#'
#' @param trt A vector of column names of treatment variables.
#' @param tau The number of time points of observation, excluding the final outcome.
#' @param time_vary A list of length tau with the column names for new time_vary to
#'  be introduced at each time point. The list should be ordered following
#'  the time ordering of the model.
#' @param baseline An optional vector of columns names for baseline covariates to be
#'  included for adjustment at every timepoint.
#' @param k An integer specifying how previous time points should be
#'  used for estimation at the given time point. Default is \code{Inf},
#'  all time points.
#'
#' @return A list of lists. Each sub-list is the same length of the
#' \code{time_vary} parameter with the variables to be used for estimation at that given time point
#' for either the treatment mechanism or outcome regression.
#' @export
#' @examples
#' \dontrun{
#' a <- c("A_1", "A_2", "A_3", "A_4")
#' bs <- c("W_1", "W_2")
#' time_vary <- list(c("L_1"), c("L_2"), c("L_3"), c("L_4"))
#'
#' # assuming no Markov property
#' create_node_list(a, 4, time_vary, bs, k = Inf)
#' }
create_node_list <- function(trt, tau, time_vary = NULL, baseline = NULL, k = Inf) {
  lifecycle::deprecate_stop("1.6.0", "create_node_list()")
}
