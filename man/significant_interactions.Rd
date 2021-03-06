% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phenopath.R
\name{significant_interactions}
\alias{significant_interactions}
\title{Significance testing for interaction features}
\usage{
significant_interactions(phenopath_fit, n = 3)
}
\arguments{
\item{phenopath_fit}{The results of a call to \code{clvm}}

\item{n}{The number of standard deviations away from 0 the posterior
estimate of beta should be to be designated significant.}
}
\value{
A logical vector describing whether each feature 
passes the significance test.
}
\description{
Given the results of \code{clvm}, decide which features show significant
iteractions between the latent trajectory and covariates. Significant 
features are designated as those where the variational mean of the interaction
coefficient falls outside the \eqn{n \sigma} interval of 0.
}
\examples{
sim <- simulate_phenopath() 
fit <- phenopath(sim$y, sim$x, elbo_tol = 1e-2)
signints <- significant_interactions(fit)
}
