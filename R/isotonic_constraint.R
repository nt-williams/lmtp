#' @importFrom OrdMonReg BoundedIsoMean
isotonic_constraint <- function(x, y) {
  na <- is.na(x) | is.na(y)
  x <- x[!na]
  y <- y[!na]

  o <- order(x)
  x <- x[o]
  y <- y[o]

  a <- rep(max(0, min(x)), length(x))
  b <- rep(max(0, min(1, max(x))), length(x))

  w <- rep(1, length(x))

  iso <- BoundedIsoMean(y, w, a, b)

  approxfun(x, iso, method = "constant", rule = 2, ties = mean)
}

cf_isotonic <- function(task, ratios, sporadic_weights, curve) {
  natural_iso <- curve$natural
  shifted_iso <- curve$shifted
  pseudo <- curve$pseudo
  for(tau in 1:task$tau) {
    for(t in 1:tau) {
      lambda <- isotonic_constraint(natural_iso[[tau]][, t], pseudo[[tau]][, t + 1])
      natural_iso[[tau]][, t] <- lambda(natural_iso[[t]][, t])
      shifted_iso[[tau]][, t] <- lambda(shifted_iso[[t]][, t])
    }
  }

  influence_functions <- matrix(0, ncol = ncol(curve$influence_functions), nrow = nrow(curve$influence_functions))
  for(t in 1:task$tau) {
    influence_functions[, t] <- eif(ratios$ratios, sporadic_weights$weights, shifted_iso[[t]], natural_iso[[t]], t = 1, tau = t)
  }

  list(
    influence_functions = influence_functions,
    shifted = shifted_iso
  )
}
