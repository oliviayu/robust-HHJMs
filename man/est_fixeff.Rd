% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/est_fixeff.r
\name{est_fixeff}
\alias{est_fixeff}
\title{Estimate fixed parameters by maximizing the profile h-likelihood function}
\usage{
est_fixeff(
  param,
  other_param,
  RespLog,
  long.data,
  surv.data,
  Bi,
  invSIGMA,
  lower = NULL,
  distribution = "normal",
  degree = 0,
  Silent = TRUE
)
}
\description{
Estimate fixed parameters by maximizing the profile h-likelihood function
}
