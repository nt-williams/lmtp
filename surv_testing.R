library(dplyr)
# library(lmtp)
library(future)
devtools::load_all()

sim_point_surv2 <- sim_point_surv
sim_point_surv2$id <- 1:nrow(sim_point_surv)
eligible_ids <- sim_point_surv2 %>% filter(Y.6 == 0) %>% pull(id)

# ivan advice, generate a time of comp event for everyone
# create indicator for competing event and modify the rest of the people

sim_point_surv3 <-
  sim_point_surv2 %>%
  mutate(CR.1 = case_when(is.na(Y.1) ~ NA_real_,
                          TRUE ~ 0),
         CR.2 = case_when(is.na(Y.2) ~ NA_real_,
                          TRUE ~ 0),
         CR.3 = case_when(id %in% eligible_ids ~ as.numeric(rbinom(nrow(.), 1, .3)),
                          is.na(Y.3) ~ NA_real_,
                          TRUE ~ 0),
         CR.4 = case_when(CR.3 == 1 ~ 1,
                          id %in% eligible_ids ~ as.numeric(rbinom(nrow(.), 1, .3)),
                          is.na(Y.4) ~ NA_real_,
                          TRUE ~ 0),
         CR.5 = case_when(CR.4 == 1 ~ 1,
                          id %in% eligible_ids ~ as.numeric(rbinom(nrow(.), 1, .3)),
                          is.na(Y.5) ~ NA_real_,
                          TRUE ~ 0),
         CR.6 = case_when(CR.5 == 1 ~ 1,
                          id %in% eligible_ids ~ as.numeric(rbinom(nrow(.), 1, .3)),
                          is.na(Y.6) ~ NA_real_,
                          TRUE ~ 0)
         )
         
  

# Example 5.1
# Time-to-event analysis with a binary time-invariant exposure. Interested in
# the effect of treatment being given to all observations on the cumulative
# incidence of the outcome.
# For a survival problem, the outcome argument now takes a vector of outcomes
# if an observation experiences the event prior to the end of follow-up, all future
# outcome nodes should be set to 1 (i.e., last observation carried forward).
A <- "trt"
Y <- paste0("Y.", 1:6)
C <- paste0("C.", 0:5)
CR <- paste0("CR.", 1:6)
W <- c("W1", "W2")
# 
#debug(lmtp_sdr)
#debug(cf_r)
# debug(cf_sdr)
#debug(at_risk)
# debug(check_at_risk)
# debug(Meta$new)

# ask Nick why I'm seeing a call for at_risk when the censoring variables are used to define at risk
# ivan thinks he may be computing risk sets for the risk ratio
# note that risk ratio can be estimated from 1-tau (not tau to 1)

progressr::with_progress(
lmtp_sdr(
  sim_point_surv3, A, Y,
   baseline = W, cens = C, comp_risk = CR, folds = 2,
  shift = static_binary_on, outcome_type = "survival"
)
)
