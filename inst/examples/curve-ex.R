\donttest{
  # Example 1.1
  # Longitudinal setting, time-varying continuous exposure bounded by 0,
  # time-varying covariates, and a binary outcome with no loss-to-follow-up.
  # Interested in the effect of a treatment policy where exposure decreases by
  # one unit at every time point if an observations observed exposure is greater
  # than or equal to 2. The true value under this intervention is about 0.305.
  head(sim_t4)

  A <- c("A_1", "A_2", "A_3", "A_4")
  L <- list(c("L_1"), c("L_2"), c("L_3"), c("L_4"))

  policy <- function(data, trt) {
    a <- data[[trt]]
    (a - 1) * (a - 1 >= 1) + a * (a - 1 < 1)
  }

  lmtp_curve(sim_t4, A, "Y", time_vary = L, shift = policy, folds = 1, mtp = TRUE)

  # Example 1.2
  # The previous example assumed that the outcome (as well as the treatment variables)
  # were directly affected by all other nodes in the past. In certain situations,
  # domain specific knowledge may suggest otherwise.
  # This can be controlled using the k argument.
  lmtp_curve(sim_t4, A, "Y", time_vary = L, shift = policy, k = 0, folds = 1, mtp = TRUE)

  # Example 1.3
  # Using the same data as examples 1.1 and 1.2, but now treating the exposure
  # as an ordered categorical variable. To account for the exposure being a
  # factor we just need to modify the shift function (and the original data)
  # so as to respect this.
  tmp <- sim_t4
  for (i in A) {
    tmp[[i]] <- factor(tmp[[i]], levels = 0:5, ordered = FALSE)
  }

  policy <- function(data, trt) {
    out <- list()
    a <- data[[trt]]
    for (i in 1:length(a)) {
      if (as.character(a[i]) %in% c("0", "1")) {
        out[[i]] <- as.character(a[i])
      } else {
        out[[i]] <- as.numeric(as.character(a[i])) - 1
      }
    }
    factor(unlist(out), levels = 0:5, ordered = FALSE)
  }

  lmtp_curve(tmp, A, "Y", time_vary = L, shift = policy, k = 0, folds = 1, mtp = TRUE)

  # Example 2.1
  # Longitudinal setting, time-varying binary treatment, time-varying covariates
  # and baseline covariates with no loss-to-follow-up. Interested in a "traditional"
  # causal effect where treatment is set to 1 at all time points for all observations.
  if (require("twang")) {
    data("iptwExWide", package = "twang")

    A <- paste0("tx", 1:3)
    W <- c("gender", "age")
    L <- list(c("use0"), c("use1"), c("use2"))

    lmtp_curve(
      iptwExWide, A, "outcome", baseline = W, time_vary = L,
      shift = static_binary_on, outcome_type = "continuous",
      mtp = FALSE, folds = 1
    )
  }

  # Example 3.1
  # Longitudinal setting, time-varying continuous treatment, time-varying covariates,
  # binary outcome with right censoring. Interested in the mean population outcome under
  # the observed exposures in a hypothetical population with no loss-to-follow-up.
  head(sim_cens)

  A <- c("A1", "A2")
  L <- list(c("L1"), c("L2"))
  C <- c("C1", "C2")
  Y <- "Y"

  lmtp_curve(sim_cens, A, Y, time_vary = L, cens = C, shift = NULL, folds = 1)

  # Example 4.1
  # Time-to-event analysis with a binary time-invariant exposure. Interested in
  # the effect of treatment being given to all observations on the cumulative
  # incidence of the outcome.
  # For a survival problem, the outcome argument now takes a vector of outcomes
  # if an observation experiences the event prior to the end of follow-up, all future
  # outcome nodes should be set to 1 (i.e., last observation carried forward).
  A <- "trt"
  Y <- paste0("Y.", 1:6)
  C <- paste0("C.", 0:5)
  W <- c("W1", "W2")

  lmtp_curve(sim_point_surv, A, Y, W, cens = C, folds = 1,
             shift = static_binary_on, outcome_type = "survival")

  # Example 5.1
  lmtp_curve(
    data = sim_competing_risks,
    trt = "A",
    cens = paste0("C", 1:5),
    compete = paste0("D", 1:5),
    baseline = paste0("W", 1:5),
    outcome = paste0("Y", 1:5),
    outcome_type = "survival",
    shift = static_binary_on,
    folds = 1
  )
}
