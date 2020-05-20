
# shift_data <- function(data, A, C, .f) {
#
#   out <- as.list(data)
#
#   if (is.null(.f)) { # only set C = 1
#     for (ce in C) {
#       out[[ce]] <- 1
#     }
#   } else {
#     for (a in A) { # shift A
#       out[[a]] <- .f(out[[a]])
#     }
#
#     for (ce in C) { # and set C = 1
#       out[[ce]] <- 1
#     }
#   }
#
#   return(as.data.frame(out))
# }

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
  return(shift_trt(shift_cens(data, cens), trt, shift))
}
