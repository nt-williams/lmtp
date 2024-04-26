cf_tmle2 <- function(task, ratios, m_init, control) {
  psi <- estimate_tmle2(task$natural,
                        task$cens,
                        task$risk,
                        task$tau,
                        m_init$mn,
                        m_init$ms,
                        ratios,
                        task$weights)

  if (!is.null(control$.boot_seed)) set.seed(control$.boot_seed)

  boots <- replicate(control$.B,
                     sample(1:nrow(task$natural), nrow(task$natural), replace = TRUE),
                     simplify = FALSE)

  Qnb <- vector("list", control$.B)
  for (i in seq_along(boots)) {
    ib <- boots[[i]]
    Qnb[[i]] <- future::future({
      psis <- estimate_tmle2(task$natural[ib, , drop = FALSE],
                             task$cens,
                             task$risk,
                             task$tau,
                             m_init$mn[ib, , drop = FALSE],
                             m_init$ms[ib, , drop = FALSE],
                             ratios[ib, , drop = FALSE],
                             task$weights[ib])
      weighted.mean(psis[, 1], task$weights[ib])
    },
    seed = TRUE)
  }

  list(psi = psi, booted = unlist(future::value(Qnb)))
}

estimate_tmle2 <- function(data, cens, risk, tau, mn, ms, ratios, weights) {
  m_eps <- matrix(nrow = nrow(data), ncol = tau + 1)
  m_eps[, tau + 1] <- data$tmp_lmtp_scaled_outcome

  fits <- vector("list", length = tau)
  for (t in tau:1) {
    i  <- censored(data, cens, t)$i
    j <- censored(data, cens, t)$j
    r <- at_risk(data, risk, t)

    wts <- ratios[i & r, t] * weights[i & r]

    fit <- sw(
      glm(
        m_eps[i & r, t + 1] ~ offset(qlogis(mn[i & r, t])),
        weights = wts,
        family = "binomial"
      )
    )

    m_eps[j & r, t] <- bound(plogis(qlogis(ms[j & r, t]) + coef(fit)))
    m_eps[!r, t] <- 0
  }

  m_eps
}
