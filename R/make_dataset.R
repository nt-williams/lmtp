make_dataset <- function(data, x, device) {
  self <- NULL
  dataset <- torch::dataset(
    name = "tmp_lmtp_dataset",
    initialize = function(data, x, device) {
      for (df in names(data)) {
        if (ncol(data[[df]]) > 0) {
          df_x <- data[[df]][, x, drop = FALSE]
          self[[df]] <- one_hot_encode(df_x) |>
            as_torch(device = device)
        }
      }
    },
    .getitem = function(i) {
      fields <- grep("data", names(self), value = TRUE)
      setNames(lapply(fields, function(x) self[[x]][i, ]), fields)
    },
    .length = function() {
      self$data$size()[1]
    }
  )
  dataset(data, x, device)
}
