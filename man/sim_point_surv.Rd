% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{sim_point_surv}
\alias{sim_point_surv}
\title{Simulated Point-treatment Survival Data}
\format{
A data frame with 2000 rows and 16 variables:
\describe{
\item{W1}{Binary baseline variable.}
\item{W2}{Categorical baseline variable.}
\item{trt}{Binary treatment variable.}
\item{C.0}{Censoring indicator that the observation is observed future time points.}
\item{Y.1}{Outcome node at time 1.}
\item{C.1}{Censoring indicator that the observation is observed future time points.}
\item{Y.2}{Outcome node at time 2.}
\item{C.2}{Censoring indicator that the observation is observed future time points.}
\item{Y.3}{Outcome node at time 3.}
\item{C.3}{Censoring indicator that the observation is observed future time points.}
\item{Y.4}{Outcome node at time 4.}
\item{C.4}{Censoring indicator that the observation is observed future time points.}
\item{Y.5}{Outcome node at time 5.}
\item{C.5}{Censoring indicator that the observation is observed future time points.}
\item{Y.6}{Final outcome node.}
}
}
\usage{
sim_point_surv
}
\description{
A dataset with a time-to-event outcome, two baseline nodes, a binary
point treatment, six past-time outcome nodes, and six censoring indicators.
}
\keyword{datasets}
