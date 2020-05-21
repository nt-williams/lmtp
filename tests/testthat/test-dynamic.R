

dynamic_vec <- function(data, trt) { # the function should either be vectorized or iterate over the rows

  if (trt == "A1") { # if the first time point set to 1
    return(rep(1, nrow(data)))
  } else { # a vectorized version
    (data[["L"]] > 0)*1 + (data[["L"]] <= 0)*0 # else return 1 or 0 based on L
  }

}

dynamic_iter <- function(data, trt) { # the function should either be vectorized or iterate over the rows
  if (trt == "A1") { # if the first time point set to 1
    return(rep(1, nrow(data)))
  } else { # an iterative version
    out <- list()
    for (i in 1:nrow(data)) {
      if (is.na(data[i, "L"])) {
        out[[i]] <- NA_real_
      } else if (data[i, "L"] > 0) {
        out[[i]] <- 1 # else return 1 or 0 based on L2
      } else {
        out[[i]] <- 0
      }
    }
    unlist(out)
  }
}

# example -----------------------------------------------------------------

library(ltmle)

rexpit <- function(x) rbinom(n=length(x), size=1, prob=plogis(x))
n <- 1000
W <- rnorm(n)
A1 <- rexpit(W)
C <- rexpit(0.6 * W - 0.5 * A1)
uncensored <- C == 1
L <- A2 <- Y <- rep(NA, n)
L[uncensored] <- (0.3 * W[uncensored] + 0.2 * A1[uncensored] + rnorm(sum(uncensored)))
A2[uncensored] <- rexpit(W[uncensored] + A1[uncensored] + L[uncensored])
Y[uncensored] <- rexpit(W[uncensored] - 0.6 * A1[uncensored] + L[uncensored] - 0.8 * A2[uncensored])
sim <- data.frame(W, A1, C, L, A2, Y)

dynamic_vec(sim, "A1")
dynamic_vec(sim, "A2")
dynamic_iter(sim, "A1")
dynamic_iter(sim, "A2")

shift_trt(sim, c("A1", "A2"), static_binary_on)
shift_trt(sim, c("A1", "A2"), static_binary_off)

shift_trt(sim, "A1", dynamic_vec)
shift_trt(sim, "A2", dynamic_vec)
shift_trt(sim, c("A1", "A2"), dynamic_vec)

shift_cens(sim, "C")
shift_data(sim, c("A1", "A2"), "C", dynamic_vec)
shift_data(sim, c("A1", "A2"), "C", NULL)
