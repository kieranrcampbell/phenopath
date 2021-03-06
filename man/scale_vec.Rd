% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clvm.R
\name{scale_vec}
\alias{scale_vec}
\title{Scale a vector to have mean 0 and variance 1}
\usage{
scale_vec(x)
}
\arguments{
\item{x}{Input vector to scale}
}
\value{
Scaled vector
}
\description{
Scales vector to mean 0 variance 1 unless input standard deviation is 0
in which case original vector is returned
}
\keyword{internal}
