% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.R
\name{simulate_phenopath}
\alias{simulate_phenopath}
\title{Simulate from phenopath model}
\usage{
simulate_phenopath(N = 100, G = 40, G_de = NULL, G_pst = NULL,
  G_pst_beta = NULL, G_de_pst_beta = NULL)
}
\arguments{
\item{N}{Number of samples to simulate}

\item{G}{Number of genes to simulate. Should be divisible by 4}

\item{G_de}{Number of genes to simulate from the \emph{differential expression} regime}

\item{G_pst}{Number of genes to simulate from the \emph{pseudotime} regime}

\item{G_pst_beta}{Number of genes to simulate from the \emph{pseudotime + beta interactions} regime}

\item{G_de_pst_beta}{Number of genes to simulate from the \emph{differential expression + pseudotime + interactions} regime}
}
\value{
A list with four entries:
\itemize{
\item \code{parameters} A \code{tibble} with an entry for each gene and a column
for each parameter value and simulation regime (see details).
\item \code{y} The N-by-G simulated gene expression matrix.
\item \code{x} The N-length vector of covariates.
\item \code{z} The N-length vector of pseudotimes.
}
}
\description{
Simulate synthetic data from the phenopath model, where each gene is sampled from 
one of four types  (see details).
}
\details{
Will simulate data for a number of genes corresponding to one of four regimes:
\enumerate{
\item \code{de} ('differential expression'), where the gene has no association to the latent 
trajectory and exhibits differential expression only
\item \code{pst} ('pseudotime'), the gene shows pseudotemporal regulation but no 
differential regulation
\item \code{pst_beta} ('pseudotime + beta interactions'), the gene shows pseudotemporal regulation
that is modulated by covariate interactions 
\item \code{de_pst_beta} ('differential expression + pseudotime + interactions), all of the above
}
}
\examples{
sim <- simulate_phenopath()
}
