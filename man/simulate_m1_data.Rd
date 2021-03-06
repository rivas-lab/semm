% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_m1_data.R
\name{simulate_m1_data}
\alias{simulate_m1_data}
\title{Simulate data for testing model 1.}
\usage{
simulate_m1_data(N, p, Sigma)
}
\arguments{
\item{N}{number of samples to generate}

\item{p}{2 d vector that contains the portion in null and then non-null components respectively}

\item{Sigma}{a 2x2 variance covariance matrix, if it is not positive definite, will use the nearest PD}
}
\value{
simulated Bs, SEs, that can be used for fitting M1
}
\description{
Simulate data for testing model 1.
}
