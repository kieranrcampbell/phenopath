


#' PhenoPath - genomic trajectories with heterogeneous backgrounds
#' 
#' PhenoPath learns genomic trajectories in the presence of heterogenous environmental
#' and genetic backgrounds. It takes input gene expression measurements that are modelled
#' by a single unobserved factor (the "trajectory"). The regulation of genes along the
#' trajectory is perturbed by an additional set of covariates (such as genetic or
#' environmental status) allowing for the identification of covariate-trajectory interactions.
#' The model is fitted using mean-field co-ordinate ascent variational inference.
#' 
#' @param exprs_obj Input gene expression, either
#' \enumerate{
#' \item An \linkS4class{ExpressionSet} object such as an \linkS4class{SCESet}, \emph{or}
#' \item A cell-by-gene matrix of normalised expression values in log form.
#' }
#' @param x The covariate vector, either
#' \enumerate{
#' \item The name of a column of \code{pData(exprs_obj)} if \code{exprs_obj} is an
#' \code{ExpressionSet}, \emph{or}
#' \item A numeric of factor vector of length equal to the number of cells
#' }
#' @param elbo_tol The relative pct change in the ELBO below which is considered converged.
#' See convergence section in details below.
#' @param z_init The initialisation of the latent trajectory. Should be one of
#' \enumerate{
#' \item A positive integer describing which principal component of the data should
#' be used for initialisation (default 1), \emph{or}
#' \item A numeric vector of length number of samples to be used directly for initialisation, \emph{or}
#' \item The text character \code{"random"}, for random initialisation from a standard
#' normal distribution.
#' }
#' @param ... Additional arguments to be passed to \code{\link{clvm}}. See description
#' below for more details or call \code{?clvm}.
#' 
#' @return An S3 structure with the following entries:
#' 
#' @details
#' \strong{Input expression}
#' 
#' If an \code{ExpressionSet} is provided, \code{exprs(...)} is used
#'
#' \strong{Encoding covariates}
#' 
#' If \code{x} is the name of...
#' 
#' \strong{Convergence of variational inference}
#' 
#' \code{plot_elbo}
#'
#' \strong{Additional arguments}
#' 
#' Addition arguments to \code{...}
#' 
#' 
#' @examples
#' sim <- simulate_phenopath() # returns a list with gene expression in y and covariates in x
#' fit <- phenopath(sim$y, sim$x, elbo_tol = 1e-2)
#' 
#' 
#' @export
phenopath <- function(exprs_obj, x, 
                      elbo_tol = 1e-5, z_init = 1, ...) {
  
  is_eset <- is_matrix <- FALSE
  N <- -1 # Number of samples
  y <- NULL # Expression object
  xx <- NULL # Covariate matrix we'll actually use
  
  # Get expression input
  if(is(exprs_obj, "ExpressionSet")) {
    is_eset <- TRUE
    N <- ncol(exprs_obj)
    y <- t(Biobase::exprs(exprs_obj))
    
  } else if(is.matrix(exprs_obj) && is.numeric(exprs_obj)) {
    is_matrix <- TRUE
    N <- nrow(exprs_obj)
    y <- exprs_obj
  } else {
    stop("Input exprs_obj must either be an ExpressionSet or numeric matrix of gene expression")
  }
  
  # Sort covariate input
  if(length(x) == 1 && is.character(x)) {
    # x describes a column of pData(exprs_obj)
    if(!is_eset || (!x %in% Biobase::varLabels(exprs_obj))) {
      stop("If x is a character then exprs_obj must be an ExpressionSet (where x represents a column in pData(exprs_obj))")
    }
    xx <- Biobase::pData(exprs_obj)[[x]]
  } else if(is.vector(x)) {
    if(length(x) != N) {
      stop("If x is a character vector it must be of compatible dimensions to exprs_obj (same number of samples)")
    }
    xx <- x
  }
  
  # Couple of extra checks
  if(is.character(xx)) {
    warning("x is a character vector...coercing to factors...")
    x <- factor(x)
  }
  if(!is.factor(xx) && !is.numeric(xx)) {
    stop("If x is a vector it must be of type numeric, factor, or character")
  }

  # If single vector then scale
  if(is.numeric(xx)) xx <- scale_vec(xx)
  
  # Now we come to the case of factors
  if(is.factor(xx)) {
    x_mat <- model.matrix(~ xx, data.frame(xx))
    x_mat <- x_mat[,-1] # Remove intercept
    x_mat <- apply(x_mat, 2, scale_vec) # Centre scale values
  } else {
    x_mat <- matrix(xx)
  }
  
  cl_fit <- clvm(y, x_mat, ...)
  
  cl_fit$N <- nrow(y)
  cl_fit$G <- ncol(y)
  
  return(
    structure(cl_fit, class = "phenopath_fit")
  )
}

#' Print a PhenoPath fit
#' 
#' @param x A \code{phenopath_fit} returned by a call to \code{phenopath}
#' @param ... Additional arguments
#' 
#' @export
#' 
#' @return A string representation of a \code{phenopath_fit} object.
#' 
#' @examples 
#' sim <- simulate_phenopath() # returns a list with gene expression in y and covariates in x
#' fit <- phenopath(sim$y, sim$x, elbo_tol = 1e-2)
#' print(fit)
print.phenopath_fit <- function(x, ...) {
  msg <- paste("PhenoPath fit with",
               x$N, "cells and", x$G, "genes")
  cat(msg)
}

#' Get the latent PhenoPath trajectory
#' 
#' @param phenopath_fit An object of class \code{phenopath_fit}
#' @return A vector of latent trajectory (pseudotime) values
#' 
#' @examples 
#' sim <- simulate_phenopath() # returns a list with gene expression in y and covariates in x
#' fit <- phenopath(sim$y, sim$x, elbo_tol = 1e-2)
#' z <- trajectory(fit)
#' 
#' @export
trajectory <- function(phenopath_fit) {
  stopifnot(is(phenopath_fit, "phenopath_fit"))
  return(phenopath_fit$m_z)
}

#' Get the interaction effect sizes
#' 
#' @param phenopath_fit An object of class \code{phenopath_fit}
#' @return TODO
#' 
#' @export
#' 
#' @examples 
#' sim <- simulate_phenopath() # returns a list with gene expression in y and covariates in x
#' fit <- phenopath(sim$y, sim$x, elbo_tol = 1e-2)
#' beta <- interaction_effects(fit)
interaction_effects <- function(phenopath_fit) {
  stopifnot(is(phenopath_fit, "phenopath_fit"))
  P <- nrow(phenopath_fit$m_beta)
  return(phenopath_fit$m_beta[1:P,,drop=TRUE])
}

#' Significance testing for interaction features
#' 
#' Given the results of \code{clvm}, decide which features show significant
#' iteractions between the latent trajectory and covariates. Significant 
#' features are designated as those where the variational mean of the interaction
#' coefficient falls outside the \eqn{n \sigma} interval of 0.
#' 
#' @param pcavi The results of a call to \code{clvm}
#' 
#' @return A logical vector describing whether each feature passes the significance test.
#' 
#' @export
significant_interactions <- function(phenopath_fit, n = 2) {
  stopifnot(is(phenopath_fit, "phenopath_fit"))
  m_beta <- phenopath_fit$m_beta
  pos_sd <- sqrt(phenopath_fit$s_beta)
  
  sig <- sapply(seq_len(nrow(m_beta)), function(i) {
    m_beta[i,] - n * pos_sd[i,] > 0 | m_beta[i,] + n * pos_sd[i,] < 0
  })
  return(sig)
}