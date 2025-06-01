# --- Helper functions for flip effect estimation --- #
define_flip_functions <- function(overlap, trimming_threshold, smoothing_constant) {
  if (overlap) {
    # Overlap weights
    flip_func       <- function(p) p * (1 - p)
    flip_func_deriv <- function(p) 1 - 2 * p
  } else {
    # Smooth trimming
    flip_func <- function(p) {
      if (trimming_threshold == 0) {
        1 - exp(-smoothing_constant * p)
      } else {
        p / (p + exp(-smoothing_constant * (p - trimming_threshold)))
      }
    }
    flip_func_deriv <- function(p) {
      if (trimming_threshold == 0) {
        smoothing_constant * exp(-smoothing_constant * p)
      } else {
        (smoothing_constant * p - 1) * exp(smoothing_constant * (trimming_threshold - p)) / 
          ((exp(trimming_threshold * smoothing_constant) + p * exp(smoothing_constant * p))^2)
      }
    }
  }
  return(list(flip_func, flip_func_deriv))
}

compute_q_phi_scores <- function(task, prop_scores, target_regime, flip_func, flip_func_deriv) {
  
  # Helper functions: set 0/0 = 0. 
  safe_divide <- function(a,b) { x <- a/b; x[is.nan(x)] <- 0; x }
  
  # Initialize q_scores and ratios
  q_scores <- matrix(NA_real_, nrow = nrow(prop_scores), ncol = ncol(prop_scores))
  ratios <- matrix(NA_real_, nrow = nrow(prop_scores), ncol = ncol(prop_scores))
  phi_scores <- array(NA_real_, dim = c(nrow(prop_scores), ncol(prop_scores), 2))
  
  for (t in 1:task$tau) {
    
    target_tx <- target_regime[t]
    A_t <- task$natural[,c(task$vars$A[[t]])]
    
    flip_p  <- flip_func(prop_scores[, t] * target_tx +
                           (1 - prop_scores[, t]) * (1 - target_tx))
    flip_p_deriv <- flip_func_deriv(prop_scores[, t] * target_tx +
                                      (1 - prop_scores[, t]) * (1 - target_tx))
    
    q_scores[, t] <- prop_scores[, t] * (1 - flip_p) + (target_tx == 1) * flip_p
    
    ratios[, t] <- 
      safe_divide(q_scores[, t] * A_t + (1 - q_scores[, t]) * (1 - A_t), 
                  prop_scores[, t] * A_t + (1 - prop_scores[, t]) * (1 - A_t))
    
    # phi_scores
    for (lvl in 0:1) {
      phi_scores[, t, lvl + 1] <- (2 * (lvl == target_tx) - 1) *
        ((A_t == target_tx) - (prop_scores[, t] * A_t + (1 - prop_scores[, t]) * (1 - A_t))) *
        (1 - flip_p + flip_p_deriv * (prop_scores[, t] * (1 - A_t) + (1 - prop_scores[, t]) * A_t))
    }
  }
  
  return(list("q_scores" = q_scores,
              "ratios" = ratios,
              "phi_scores" = phi_scores))
}