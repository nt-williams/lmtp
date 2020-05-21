#
#
# dynamic_vec <- function(data, trt) { # the function should either be vectorized or iterate over the rows
#
#   if (trt == "A1") { # if the first time point set to 1
#     return(rep(1, nrow(data)))
#   } else { # a vectorized version
#     (data[["L"]] > 0)*1 + (data[["L"]] <= 0)*0 # else return 1 or 0 based on L
#   }
#
# }
#
# dynamic_iter <- function(data, trt) { # the function should either be vectorized or iterate over the rows
#   if (trt == "A1") { # if the first time point set to 1
#     return(rep(1, nrow(data)))
#   } else { # an iterative version
#     out <- list()
#     for (i in 1:nrow(data)) {
#       if (is.na(data[i, "L"])) {
#         out[[i]] <- NA_real_
#       } else if (data[i, "L"] > 0) {
#         out[[i]] <- 1 # else return 1 or 0 based on L2
#       } else {
#         out[[i]] <- 0
#       }
#     }
#     unlist(out)
#   }
# }
#
# # example -----------------------------------------------------------------
#
# library(ltmle)
#
# rexpit <- function(x) rbinom(n=length(x), size=1, prob=plogis(x))
# n <- 10000
# W <- rnorm(n)
# A1 <- rexpit(W)
# C1 <- rexpit(0.6 * W - 0.5 * A1)
# uncensored <- C1 == 1
# L <- A2 <- C2 <- Y <- rep(NA, n)
# L[uncensored] <- (0.3 * W[uncensored] + 0.2 * A1[uncensored] + rnorm(sum(uncensored)))
# A2[uncensored] <- rexpit(W[uncensored] + A1[uncensored] + L[uncensored])
# C2[uncensored] <- 1
# C2[!uncensored] <- 0
# Y[uncensored] <- rexpit(W[uncensored] - 0.6 * A1[uncensored] + L[uncensored] - 0.8 * A2[uncensored])
# sim <- data.frame(W, A1, C1, L, A2, C2, Y)
#
# dynamic_vec(sim, "A1")
# dynamic_vec(sim, "A2")
# dynamic_iter(sim, "A1")
# dynamic_iter(sim, "A2")
#
# shift_trt(sim, c("A1", "A2"), static_binary_on)
# shift_trt(sim, c("A1", "A2"), static_binary_off)
#
# shift_trt(sim, "A1", dynamic_vec)
# shift_trt(sim, "A2", dynamic_vec)
# shift_trt(sim, c("A1", "A2"), dynamic_vec)
#
# shift_cens(sim, "C")
# shift_data(sim, c("A1", "A2"), "C", dynamic_vec)
# shift_data(sim, c("A1", "A2"), "C", NULL)
#
# # data generation
# a <- c("A1", "A2")
# baseline <- "W"
# cens <- c("C1", "C2")
# nodes <- list(c(NULL), c("L"))
# lrnrs <- sl3::make_learner_stack(sl3::Lrnr_glm,
#                                  sl3::Lrnr_mean,
#                                  sl3::Lrnr_xgboost)
#
# # estimators
# # static on off at all time points
# progressr::with_progress({
#   psi1 <- lmtp_tmle(sim, a, "Y", nodes, baseline, cens, shift = static_binary_on,
#             outcome_type = "binomial", folds = 2,
#             learners_outcome = lrnrs, learners_trt = lrnrs)
# })
#
# lmtp_tmle(sim, a, "Y", nodes, baseline, cens, shift = static_binary_on,
#             outcome_type = "binomial", folds = 2,
#           learners_outcome = lrnrs, learners_trt = lrnrs)
#
# lmtp_tmle(sim, a, "Y", nodes, baseline, cens, shift = static_binary_off,
#           outcome_type = "binomial", folds = 2)
#
# # on at 1 off at 2
# time_vary_on <- function(data, trt) {
#   if (trt == "A1") return(rep(1, nrow(data)))
#   else return(rep(0, nrow(data)))
# }
#
# # truth = 0.433
# progressr::with_progress({
#   psi.tv <- lmtp_tmle(sim, a, "Y", nodes, baseline, cens, shift = time_vary_on,
#                     outcome_type = "binomial", folds = 5,
#                     learners_outcome = lrnrs, learners_trt = lrnrs)
# })
#
# progressr::with_progress({
#   sdr.tv <- lmtp_sdr(sim, a, "Y", nodes, baseline, cens, shift = time_vary_on,
#                       outcome_type = "binomial", folds = 5,
#                       learners_outcome = lrnrs, learners_trt = lrnrs)
# })
#
# # truth = 0.345
# progressr::with_progress({
#   psi.fme <- lmtp_tmle(sim, a, "Y", nodes, baseline, cens, shift = dynamic_vec,
#                       outcome_type = "binomial", folds = 5,
#                       learners_outcome = lrnrs, learners_trt = lrnrs)
# })
#
# progressr::with_progress({
#   sdr.fme <- lmtp_sdr(sim, a, "Y", nodes, baseline, cens, shift = dynamic_vec,
#                        outcome_type = "binomial", folds = 5,
#                        learners_outcome = lrnrs, learners_trt = lrnrs)
# })
#