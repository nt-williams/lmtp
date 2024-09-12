#' Sequential neural network module function factory
#'
#' @param layers \[numeric(1)\]\cr Number of hidden layers.
#' @param hidden \[numeric(1)\]\cr Number of hidden units.
#' @param dropout \[numeric(1)\]\cr Dropout rate.
#'
#' @return A function that returns a sequential neural network module.
#' @export
#'
#' @examples
#' if (torch::torch_is_installed()) sequential_module()
sequential_module <- function(layers = 1, hidden = 20, dropout = 0.1) {
  function(d_in) {
    d_out <- 1

    middle_layers <- lapply(1:layers, \(x) torch::nn_sequential(torch::nn_linear(hidden, hidden), torch::nn_elu()))

    torch::nn_sequential(
      torch::nn_linear(d_in, hidden),
      torch::nn_elu(),
      do.call(torch::nn_sequential, middle_layers),
      torch::nn_linear(hidden, d_out),
      torch::nn_dropout(dropout),
      torch::nn_softplus()
    )
  }
}
