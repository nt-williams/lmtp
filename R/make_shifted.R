make_shifted <- function(data, A, C, shift, shifted) {
  if (!is.null(shifted)) {
    return(shifted)
  }

  if (is.null(shifted) && !is.null(shift)) {
    return(shift_data(data, A, C, shift))
  }

  if (is.null(shifted) && is.null(shift)) {
    return(shift_data(data, A, C, shift))
  }
}
