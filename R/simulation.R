
#' Simulate one gene
#' 
#' Simulate one gene from the model given paramters, z and covariates
#' 
#' @keywords internal
simulate_one_gene <- function(N, pst, x, alpha = 0, c = 0, beta = 0, tau = 1e6) {
  mu <- alpha * x + (c + beta * x) * pst
  rnorm(N, mu, 1 / sqrt(tau))
}

#' Sample parameters for simulation
#' 
#' Sample parameters from de regime
#' @keywords internal
#' @name sample_fns
sample_de <- function() {
  alpha <- sample(c(-1, 1), 1)
  c <- beta <- 0
  return(c(alpha, c, beta))
}

#' Sample parameters from pseudotime regime
#' @keywords internal
#' @name sample_fns
sample_pst <- function() {
  c <- 1 * sample(c(-1, 1), 1)
  # c <- rnorm(1)
  alpha <- beta <- 0
  return(c(alpha, c, beta))
}

#' Sample parameters from pseudotime-interaction regime
#' @keywords internal
#' @name sample_fns
sample_pst_beta <- function() {
  c <- 1 * sample(c(-1, 1), 1)  
  # c <- rnorm(1)
  beta <- 1 * sample(c(-1, 1), 1)  
  alpha <- 0
  return(c(alpha, c, beta))
}

#' Sample paramters from de-pseudotime-interaction regime
#' @keywords internal
#' @name sample_fns
sample_de_pst_beta <- function() {
  c <- 1 * sample(c(-1, 1), 1)  
  # c <- rnorm(1)
  beta <- 1 * sample(c(-1, 1), 1)  
  alpha <- sample(c(-1, 1), 1)
  return(c(alpha, c, beta))
}


#' Simulate from phenopath model
#' 
#' Simulate synthetic data from the phenopath model, where each gene is sampled from 
#' one of four types  (see details).
#' 
#' @param N Number of samples to simulate
#' @param G Number of genes to simulate. Should be divisible by 4
#' @param G_de Number of genes to simulate from the \emph{differential expression} regime
#' @param G_pst Number of genes to simulate from the \emph{pseudotime} regime
#' @param G_pst_beta Number of genes to simulate from the \emph{pseudotime + beta interactions} regime
#' @param G_de_pst_beta Number of genes to simulate from the \emph{differential expression + pseudotime + interactions} regime
#' 
#' @return A list with four entries:
#' \itemize{
#' \item \code{parameters} A \code{tibble} with an entry for each gene and a column
#' for each parameter value and simulation regime (see details).
#' \item \code{y} The N-by-G simulated gene expression matrix.
#' \item \code{x} The N-length vector of covariates.
#' \item \code{z} The N-length vector of pseudotimes.
#' }
#' 
#' 
#' @details
#' 
#' Will simulate data for a number of genes corresponding to one of four regimes:
#' \enumerate{
#' \item \code{de} ('differential expression'), where the gene has no association to the latent 
#' trajectory and exhibits differential expression only
#' \item \code{pst} ('pseudotime'), the gene shows pseudotemporal regulation but no 
#' differential regulation
#' \item \code{pst_beta} ('pseudotime + beta interactions'), the gene shows pseudotemporal regulation
#' that is modulated by covariate interactions 
#' \item \code{de_pst_beta} ('differential expression + pseudotime + interactions), all of the above
#' }
#' 
#' @export
#' @importFrom dplyr mutate
#' @importFrom tibble as_data_frame
#' 
#' @examples
#' sim <- simulate_phenopath()
simulate_phenopath <- function(N = 100, G = 40,
                               G_de = NULL, G_pst = NULL,
                               G_pst_beta = NULL, G_de_pst_beta = NULL) {

  setifnotnull <- function(x) ifelse(is.null(x), G / 4, x)
  if(any(is.null(c(G_de, G_pst, G_pst_beta, G_de_pst_beta))) && (G %% 4 != 0)) {
    stop("If any of G_de, G_pst, G_pst_beta, G_de_pst_beta are left NULL then G must be divisible by 4 four default")
  }
  if(all(!is.null(c(G_de, G_pst, G_pst_beta, G_de_pst_beta)))) {
    G <- G_de + G_pst + G_pst_beta + G_de_pst_beta
  }
  G_de <- setifnotnull(G_de)
  G_pst <- setifnotnull(G_pst)
  G_pst_beta <- setifnotnull(G_pst_beta)
  G_de_pst_beta <- setifnotnull(G_de_pst_beta)
  
  
  x <- rep(c(1, -1), each = N / 2)
  pst <- rnorm(N)
  pst <- (pst - mean(pst)) / sd(pst)

  de_pars <- pst_pars <- pst_beta_pars <- de_pst_beta_pars <- NULL
  if(G_de > 0)
    de_pars <- t(replicate(G_de, sample_de()))
  if(G_pst > 0)
    pst_pars <- t(replicate(G_pst, sample_pst()))
  if(G_pst_beta > 0)
    pst_beta_pars <- t(replicate(G_pst_beta, sample_pst_beta()))
  if(G_de_pst_beta > 0)
    de_pst_beta_pars <- t(replicate(G_de_pst_beta, sample_de_pst_beta()))
  
  regimes <- c("de", "pst", "pst_beta", "de_pst_beta")
  
  par_list <- list(de_pars, pst_pars, pst_beta_pars, de_pst_beta_pars)
  par_list <- par_list[!sapply(par_list, is.null)]
  pars <- as_data_frame(do.call("rbind", par_list))
  names(pars) <- c("alpha", "lambda", "beta")
  pars <- mutate(pars, regime = rep(regimes, times = c(G_de, G_pst, G_pst_beta, G_de_pst_beta)))
  
  #pars$alpha <- scale_vec(pars$alpha)
  #pars$lambda <- scale_vec(pars$lambda)
  #pars$beta <- scale_vec(pars$beta)
  
  gex <- apply(pars, 1, function(p) {
    p <- as.numeric(p[1:3])
    simulate_one_gene(N, pst, x, p[1], p[2], p[3])
  })
  
  return(list(parameters = pars, y = gex, x = x, z = pst))
}