% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sensitivity_functions.R
\name{make_prior}
\alias{make_prior}
\title{Construct prior list based on type of sensitivity analysis}
\usage{
make_prior(params, prior_type = c("unit", "partial", "strict", "dd",
  "powerprior"), dd_control = NULL)
}
\arguments{
\item{params}{List of parameters}

\item{prior_type}{Variance type: unit variance ("unit"), informative priors 
on non-identifiable parameters ("partial"), strong prior information 
("strict"), data-driven ("dd") based on external or simulated data set, or
data pooling/power prior ("pp") based on the external or simulated data set}

\item{dd_control}{Named list of data-driven options. If external data is to
be used, "small_data" needs to be formatted as output from 
\code{\link{simulate_data}}. If secondary data source is to be simulated, 
need a list of parameters as from \code{\link{return_dgp_parameters}}.}
}
\value{
Named list of prior means and variance-covariance matrices
}
\description{
Construct prior list based on type of sensitivity analysis
}
