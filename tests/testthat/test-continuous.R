

library(lmtp)
library(sl3)
library(future)
library(progressr)

pta <- readRDS("/Volumes/Research_share/Nick\ Williams/ctsc/alexiades_ctsc/alexiades_hearing_aids/use-data.rds")

a <- paste0("aided_yet_", 0:4)
nodes <- list(c("init_age", "gender", "base_value"),
              c("value_0"),
              c("value_1"),
              c("value_2"),
              c("value_3"))
cens <- paste0("cens_", 0:4)
d <- function(x) 1
lrnrs <- make_learner_stack(Lrnr_mean,
                            Lrnr_glm,
                            Lrnr_earth_caret,
                            Lrnr_xgb_caret)

plan(multiprocess)
with_progress({
  fit <-
    lmtp_tmle(pta, a, "value_4", nodes, cens = cens,
              outcome_type = "continuous", shift = d,
              learners_outcome = lrnrs,
              learners_trt = lrnrs, folds = 10)
})

with_progress({
  fit0 <-
    lmtp_tmle(pta, a, "value_4", nodes, cens = cens,
              outcome_type = "continuous", shift = shift_0,
              learners_outcome = lrnrs,
              learners_trt = lrnrs, folds = 10)
})

with_progress({
  sdr <-
    lmtp_sdr(pta, a, "value_4", nodes, cens = cens,
             outcome_type = "continuous", shift = d,
             learners_outcome = lrnrs,
             learners_trt = lrnrs, folds = 10)
})

with_progress({
  sdr0 <-
    lmtp_sdr(pta, a, "value_4", nodes, cens = cens,
             outcome_type = "continuous", shift = shift_0,
             learners_outcome = lrnrs,
             learners_trt = lrnrs, folds = 10)
})


# set exposure to 1 for all observations at all time points
shift <- function(x) 1

# set exposure to 0 for all observations at all time points
shift_0 <- function(x) 0

# a binary trt example
data("iptwExWide", package = "twang")

a <- paste0("tx", 1:3)
nodes <- list(c("gender", "age", "use0"),
              c("use1"),
              c("use2"))

sdr1 <-
  lmtp_sdr(iptwExWide, a, "outcome", nodes,
           outcome_type = "continuous", shift = shift,
           learners_outcome = lrnrs,
           learners_trt = lrnrs, folds = 10)

sdr0 <-
  lmtp_sdr(iptwExWide, a, "outcome", nodes,
           outcome_type = "continuous", shift = shift_0,
           learners_outcome = lrnrs,
           learners_trt = lrnrs, folds = 10)


tmle1 <-
  lmtp_tmle(iptwExWide, a, "outcome", nodes,
           outcome_type = "continuous", shift = shift,
           learners_outcome = lrnrs,
           learners_trt = lrnrs, folds = 10)
lrnrs <- make_learner_stack(Lrnr_mean, Lrnr_glm, Lrnr_earth_caret)

lmtp_tmle(iptwExWide, a, "outcome", nodes, shift = shift,
          outcome_type = "continuous", folds = 10)
