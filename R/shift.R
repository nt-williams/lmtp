
shift_cens <- function(data, cens) {
  out <- as.list(data)
  for (ce in cens) {
    out[[ce]] <- 1
  }
  return(as.data.frame(out))
}

shift_trt <- function(data, trt, .f) {
  for (a in trt) {
    data[[a]] <- .f(data, a)
  }
  return(data)
}

shift_data <- function(data, trt, cens, shift) {
  # NEED TO FIGURE OUT WHAT TO DO IF SHIFT IS NULL
  return(shift_trt(shift_cens(data, cens), trt, shift))
}
