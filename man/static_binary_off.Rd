% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shift.R
\name{static_binary_off}
\alias{static_binary_off}
\title{Turn All Treatment Nodes Off}
\usage{
static_binary_off(data, trt)
}
\arguments{
\item{data}{A dataframe containing the treatment variables.}

\item{trt}{The name of the current treatment variable.}
}
\value{
A dataframe with all treatment nodes set to 0.
}
\description{
A pre-packaged shift function for use with provided estimators when the exposure is binary.
Used to estimate the population intervention effect when all treatment variables are set to 0.
}
\examples{
\donttest{
data("iptwExWide", package = "twang")
a <- paste0("tx", 1:3)
baseline <- c("gender", "age")
tv <- list(c("use0"), c("use1"), c("use2"))
lmtp_sdr(iptwExWide, a, "outcome", baseline = baseline, time_vary = tv,
         shift = static_binary_off, outcome_type = "continuous", folds = 2)
}
}
\seealso{
\code{\link[=lmtp_tmle]{lmtp_tmle()}}, \code{\link[=lmtp_sdr]{lmtp_sdr()}}
}
