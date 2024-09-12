cf_G <- function(task, learners, mtp, control, pb) {
  out <- vector("list", length = length(task$folds))
  for (fold in seq_along(task$folds)) {
    out[[fold]] <- future::future({
      estimate_G(
        get_folded_data(task$natural, task$folds, fold),
        get_folded_data(task$shifted, task$folds, fold),
        get_folded_data(task$conditional, task$folds, fold),
        task$trt,
        task$cens,
        task$risk,
        task$tau,
        task$node_list$trt,
        learners,
        pb,
        mtp,
        control
      )
    },
    seed = TRUE)
  }

  out <- future::value(out)

  list(
    G = recombine_outcome(out, "G", task$folds),
    fits = lapply(out, function(x) x[["fits"]])
  )
}

estimate_G <- function() {

}
