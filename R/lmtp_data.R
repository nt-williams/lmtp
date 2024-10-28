as_lmtp_wide_data <- function(data, A, Y, W, L, C, id) {
  assert_lmtp_data(data, A, Y, W, L, C, id)

  data <- data.table::copy(as.data.frame(data))
  data$lmtp_id <- create_ids(data, id)
  data <- fix_censoring_ind(data, C)
}

as_lmtp_long_data <- function() {

}
