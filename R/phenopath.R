


#' PhenoPath - genomic trajectories with heterogeneous backgrounds
#' 
#' PhenoPath learns genomic trajectories in the presence of heterogenous environmental
#' and genetic backgrounds using scalable Bayesian variational inference.
#' 
#' @param exprs_obj Input gene expression, either
#' \enumerate{
#' \item An \code{ExpressionSet} object such as an \code{SCESet}, \emph{or}
#' \item A cell-by-gene matrix of normalised expression values in log form.
#' }
#' @param x The covariate vector, either
#' \enumerate{
#' \item The name of a column of \code{pData(exprs_obj)} if \code{exprs_obj} is an
#' \code{ExpressionSet}, \emph{or}
#' \item A numeric of factor vector of length equal to the number of cells
#' }
#' 
#' 
#' @details
#' 
#' If an \code{ExpressionSet} is provided, \code{exprs(...)} is used
#' 
#' If \code{x} is the name of...
#' 
#' 
#' @export
phenopath <- function(exprs_obj, x, ...) {
  
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
  
  clvm(y, x_mat, ...)
}