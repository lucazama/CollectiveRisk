% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Approx.R
\name{simulator}
\alias{simulator}
\title{Simulator of random numbers from some not-so-common distributions}
\usage{
simulator(dist, par, n)
}
\arguments{
\item{dist}{distribution from which to simulate c("Pareto", "TruncWeibull", TruncLogNormal")}

\item{par}{parameters of the distribution chosen.}

\item{n}{the number of simulations.}
}
\value{
The function returns the simulations of Y. If you want to compute the probability P(Y<m)
you can simply run the following command mean(Y<m)
}
\description{
The function simulates random numbers from some distributions
}
\note{
Confidence intervals still to be added (along with more distributions)
}
