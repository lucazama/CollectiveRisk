% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Moments.R
\name{rawtostandard}
\alias{rawtostandard}
\title{Function to convert raw Moments to Standard ones.}
\usage{
rawtostandard(moments)
}
\arguments{
\item{moments}{is a vector of raw moments.}
}
\value{
A vector of standardized moments. the nth element of the vector
is the nth standardized moment.
}
\description{
The function converts the raw moments of a distribution to the standardized ones.
}
\note{
Maximum length supported is n=5.
}
