eif_mtp <- function(r, shifted, natural, t, tau) {
  if (missing(tau)) tau <- ncol(r)
  if (missing(t)) t <- 1
  # natural[is.na(natural)] <- -999
  # shifted[is.na(shifted)] <- -999
  m <- shifted[, (t + 1):(tau + 1), drop = FALSE] - natural[, t:tau, drop = FALSE]
  rowSums(compute_weights(r, t, tau) * m, na.rm = TRUE) + shifted[, t]
}

compute_weights <- function(r, t, tau) {
  out <- t(apply(r[, t:tau, drop = FALSE], 1, cumprod))
  if (ncol(out) > ncol(r)) return(t(out))
  out
}

#' Internal: Calculates efficient influence function with estimated intervention
#' propensity scores
eif_Q <- function(task, data, seq_ests, q_scores, phi_scores, ratios, time, final_time) {
  
  # plug-in estimate of pseudo-outcome
  pseudo <- seq_ests[, time, 1] * (1 - q_scores[, time]) + 
    seq_ests[, time, 2] * q_scores[, time]

  pseudo <- pseudo + rowSums(seq_ests[, time, ] * phi_scores[, time, ])
    
  for (s in time:final_time) {
    
    for (k in time:s) {
      if (k == time) {
        ratio_product <- ratios[, k]
        
      } else {
        ratio_product <- ratio_product * ratios[, k]
      }
    }
    
    if (s < final_time) {
      resid <- 
        seq_ests[, s+1, 1] * (1 - q_scores[, s+1]) + seq_ests[, s+1, 2] * q_scores[, s+1] +
        rowSums(seq_ests[, s+1, ] * phi_scores[, s+1, ]) -
        seq_ests[, s, ][cbind(1:nrow(data), data[, task$vars$A[[time]]] + 1)]
    } else { 
      # Use final outcome data for final residual
      resid <- data[, c("final")] -
        seq_ests[, s, ][cbind(1:nrow(data), data[, task$vars$A[[time]]]+1)]
    }
    
    pseudo <- pseudo + ratio_product * resid
  }
  
  return(pseudo)
}
