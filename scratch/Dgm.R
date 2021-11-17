R6::R6Class(
  "Dgm",
  public = list(
    initialize = function(n) {

    }
  ),
  private = list(
    L = function(a1, a2, l) {
      plogis(-0.3 * l + 0.5 * a1 - 0.1 * a2)
    },
    A1 = function(prev_a1, prev_a2, l) {

    },
    A2 = function(prev_a1, prev_a2, l) {

    },
    Y = function(a1, a2, l) {
      plogis(-2 + 1 / (1 - 1.2 * a1 - 0.3 * a2 - 0.1 * l))
    }
  )
)
