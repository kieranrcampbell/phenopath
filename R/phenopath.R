


#' PhenoPath - genomic trajectories with heterogeneous backgrounds
#' 
#' PhenoPath learns genomic trajectories in the presence of 
#' heterogenous environmental and genetic backgrounds. It takes 
#' input gene expression measurements that are modelled
#' by a single unobserved factor (the "trajectory"). The regulation 
#' of genes along the trajectory is perturbed by an additional set of 
#' covariates (such as genetic or environmental status) allowing for 
#' the identification of covariate-trajectory interactions.
#' The model is fitted using mean-field co-ordinate ascent variational inference.
#' 
#' @param exprs_obj Input gene expression, either
#' \enumerate{
#' \item An \linkS4class{ExpressionSet} object, \emph{or}
#' \item A cell-by-gene matrix of normalised expression values in log form.
#' }
#' @param x The covariate vector, either
#' \enumerate{
#' \item The name of a column of \code{pData(exprs_obj)} if \code{exprs_obj} is an
#' \code{ExpressionSet}, \emph{or}
#' \item A numeric of factor vector of length equal to the 
#' number of cells, \emph{or}
#' \item A formula from which to build a model matrix from \code{pData(exprs_obj)}, 
#' if \code{exprs_obj} is a \linkS4class{ExpressionSet}
#' }
#' @param elbo_tol The relative pct change in the ELBO below 
#' which is considered converged.
#' See convergence section in details below.
#' @param z_init The initialisation of the latent trajectory. Should be one of
#' \enumerate{
#' \item A positive integer describing which principal component of the data should
#' be used for initialisation (default 1), \emph{or}
#' \item A numeric vector of length number of samples to be used 
#' directly for initialisation, \emph{or}
#' \item The text character \code{"random"}, for random initialisation 
#' from a standard normal distribution.
#' }
#' @param ... Additional arguments to be passed to \code{\link{clvm}}. 
#' See description below for more details or call \code{?clvm}.
#' 
#' @return An S3 structure with the following entries:
#' \itemize{
#' \item \code{m_z} The converged mean estimates of the trajectory
#' \item \code{s_z} The converged standard deviation estimates of z
#' \item \code{m_beta} A P-by-G matrix of interaction coefficients
#' \item \code{s_beta} A P-by-G matrix of interaction standard deviations
#' }
#' 
#' @details
#' \strong{Input expression}
#' 
#' If an \code{ExpressionSet} is provided, \code{exprs(...)} is used. 
#' This is assumed to be in
#' a form that is suitably normalised and approximately normal, such as 
#' the log of TPM values (plus a suitable offset) or 
#' similar.
#'
#' \strong{Encoding covariates}
#' 
#' If \code{x} is the name of...
#' 
#' \strong{Convergence of variational inference}
#' 
#' It is strongly recommended to call \code{plot_elbo(phenopath_fit)} after the fitting procedure to
#' ensure the ELBO has approximately converged (though convergence metrics are printed during the
#' fitting process). For a good introduction to variational inference see Blei, D.M., Kucukelbir, A. & McAuliffe, J.D., 2017. Variational Inference: A Review for Statisticians. Journal of the American Statistical Association.
#'
#' \strong{Additional arguments}
#' 
#' Addition arguments to \code{clvm} are passed via \code{...}. For full documentation, call \code{?clvm}.
#' Some notable options:
#' \itemize{
#' \item \code{thin} - The ELBO is expensive to compute for larger datasets. The model will
#' compute the ELBO and compare convergence every \code{thin} iterations.
#' \item \code{q} and \code{tau_q} - Priors (such as capture times) for the latent space. Note that
#' \code{model_mu} should be true if \code{q} is non-zero.
#' \item \code{scale_y} By default the input expression is centre-scaled for each gene. If \code{scale_y}
#' is \code{FALSE} this does not happen - but note that \code{model_mu} should be \code{TRUE} in such a case.
#' }
#' 
#' @seealso \code{\link{clvm}} for the underlying CAVI function, \code{\link{trajectory}}
#' to extract the latent trajectory, \code{\link{interaction_effects}} for the interaction effect
#' sizes, \code{\link{significant_interactions}} for the results of Bayesian significance testing.
#' 
#' @examples
#' sim <- simulate_phenopath() # returns a list with gene expression in y and covariates in x
#' fit <- phenopath(sim$y, sim$x, elbo_tol = 1e-2)
#' 
#' # Extract the trajectory
#' z <- trajectory(fit)
#' 
#' @importFrom methods is
#' @importFrom stats model.matrix
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
  } else if(is.vector(x) | is.factor(x)) {
    if(length(x) != N) {
      stop("If x is a character vector or factor it must be of compatible dimensions to exprs_obj (same number of samples)")
    }
    xx <- x
  } else if(is(x, "formula")) {
    xx <- x
  } else if(is.matrix(x)) {
    stopifnot(nrow(x) == N)
    x_mat <- xx <- x
  }
  
  # Couple of extra checks
  if(is.character(xx)) {
    warning("x is a character vector...coercing to factors...")
    xx <- factor(xx)
  }
  if(!is.factor(xx) && !is.numeric(xx) && !is(xx, "formula")) {
    stop("If x is a vector it must be of type numeric, factor, or character")
  }

  # If single vector then scale
  if(is.numeric(xx) && !is.matrix(xx)) xx <- scale_vec(xx)
  
  # Now we come to the case of factors
  if(is.factor(xx)) {
    x_mat <- model.matrix(~ xx, data.frame(xx))
    # COME BACK AND CHANGE ME
    x_mat <- x_mat[,-1, drop = FALSE] # Remove intercept
    x_mat <- apply(x_mat, 2, scale_vec) # Centre scale values
    # x_mat <- apply(x_mat, 2, function(x) return(2*x - 1))
  } else if(is(xx, "formula")) {
    if(!is(exprs_obj, "ExpressionSet")) {
      stop("If x is a formula, y must be an ExpressionSet")
    }
    x_mat <- model.matrix(xx, Biobase::pData(exprs_obj))
    x_mat <- x_mat[,-1, drop = FALSE] # Remove intercept
    x_mat <- apply(x_mat, 2, scale_vec) # Centre scale values
  } else if(is.vector(xx)) {
    x_mat <- matrix(xx)
  }
  
  cl_fit <- clvm(y, x_mat, elbo_tol = elbo_tol, z_init = z_init, ...)
  
  cl_fit$N <- nrow(y)
  cl_fit$G <- ncol(y)
  
  feature_names <- colnames(y)
  if(is.null(feature_names)) feature_names <- paste0("feature_", seq_len(cl_fit$G))
  cl_fit$feature_names <- feature_names
  
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
#' sim <- simulate_phenopath() 
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
#' sim <- simulate_phenopath() 
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
#' sim <- simulate_phenopath() 
#' fit <- phenopath(sim$y, sim$x, elbo_tol = 1e-2)
#' beta <- interaction_effects(fit)
interaction_effects <- function(phenopath_fit) {
  stopifnot(is(phenopath_fit, "phenopath_fit"))
  P <- nrow(phenopath_fit$m_beta)
  return(phenopath_fit$m_beta[1:P,,drop=TRUE])
}

#' Get the interaction standard deviations
#' 
#' @param phenopath_fit An object of class \code{phenopath_fit}
#' @return TODO
#' 
#' @export
#' 
#' @examples 
#' sim <- simulate_phenopath() 
#' fit <- phenopath(sim$y, sim$x, elbo_tol = 1e-2)
#' beta_sd <- interaction_sds(fit)
interaction_sds <- function(phenopath_fit) {
  stopifnot(is(phenopath_fit, "phenopath_fit"))
  P <- nrow(phenopath_fit$m_beta)
  return(sqrt(phenopath_fit$s_beta[1:P,,drop=TRUE]))
}

#' Significance testing for interaction features
#' 
#' Given the results of \code{clvm}, decide which features show significant
#' iteractions between the latent trajectory and covariates. Significant 
#' features are designated as those where the variational mean of the interaction
#' coefficient falls outside the \eqn{n \sigma} interval of 0.
#' 
#' @param phenopath_fit The results of a call to \code{clvm}
#' @param n The number of standard deviations away from 0 the posterior
#' estimate of beta should be to be designated significant.
#' 
#' @return A logical vector describing whether each feature 
#' passes the significance test.
#' 
#' @export
#' @examples 
#' sim <- simulate_phenopath() 
#' fit <- phenopath(sim$y, sim$x, elbo_tol = 1e-2)
#' signints <- significant_interactions(fit)
significant_interactions <- function(phenopath_fit, n = 3) {
  stopifnot(is(phenopath_fit, "phenopath_fit"))
  m_beta <- phenopath_fit$m_beta
  pos_sd <- sqrt(phenopath_fit$s_beta)
  
  sig <- sapply(seq_len(nrow(m_beta)), function(i) {
    m_beta[i,] - n * pos_sd[i,] > 0 | m_beta[i,] + n * pos_sd[i,] < 0
  })
  if(ncol(sig) == 1) sig <- as.vector(sig)
  return(sig)
}


#' Tidy summary of interactions
#' 
#' Computes a tidy data frame of the interaction effects, pathway scores, and significance for
#' each gene and covariate interaction.
#' 
#' @param phenopath_fit An object returned by a call to \code{phenopath}
#' @param n The number of standard deviations away from 0 the posterior of the 
#' interaction effect sizes should be to be designated as significant
#' 
#' @importFrom tidyr gather
#' @export
#' 
#' @return A data frame with the following columns:
#' \itemize{
#' \item \code{feature} The feature (usually gene)
#' \item \code{covariate} The covariate, specified from the order originally supplied to 
#' the call to \code{phenopath}
#' \item \code{interaction_effect_size} The effect size of the interaction (\eqn{beta} 
#' from the statistical model)
#' \item \code{significant} Boolean for whether the interaction effect is significantly
#' different from 0
#' \item \code{chi} The precision of the ARD prior on \eqn{beta}
#' \item \code{pathway_score} The pathway score \eqn{lambda}, showing the overall
#' effect for each gene marginalised over the covariate
#' }
#' @examples 
#' sim <- simulate_phenopath() 
#' fit <- phenopath(sim$y, sim$x, elbo_tol = 1e-2)
#' interactions(fit)
interactions <- function(phenopath_fit, n = 3) {
  covariate <- interaction_effect_size <- feature <- significant_interaction <- NULL
  
  ie <- interaction_effects(phenopath_fit)
  sig <- significant_interactions(phenopath_fit)
  if(is.vector(ie)) {
    ie <- matrix(ie)
  }
  if(is.vector(sig)) {
    sig <- matrix(sig)
  }
  
  chi <- dplyr::as_data_frame(t(phenopath_fit$chi_exp))
  
  ie <- dplyr::as_data_frame(ie)
  sig <- dplyr::as_data_frame(sig)
  names(ie) <- names(sig) <- names(chi) <- paste0("covariate_", 1:ncol(ie))
  ie$feature <- sig$feature <- chi$feature <- phenopath_fit$feature_names
  
  ie_tidy <- gather(ie, covariate, interaction_effect_size, -feature)  
  sig_tidy <- gather(sig, covariate, significant_interaction, -feature)  
  chi_tidy <- gather(chi, covariate, chi, -feature)
  
  interaction_df <- dplyr::inner_join(ie_tidy, sig_tidy, by = c("feature", "covariate"))
  interaction_df <- dplyr::inner_join(interaction_df, chi_tidy, by = c("feature", "covariate"))
  
  lambda_df <- dplyr::data_frame(pathway_score = phenopath_fit$m_lambda,
                          feature = phenopath_fit$feature_names)
  
  interaction_df <- dplyr::inner_join(interaction_df, lambda_df, by = "feature")
  interaction_df
}




