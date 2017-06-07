library(Rcpp)


#' Fit a CLVM Model
#' 
#' Fit a covariate latent variable model using coordinate ascent
#' variational inference.
#' 
#' @param y A N-by-G (dynamic) input matrix 
#' @param x A N-by-P (static) input matrix
#' @param maxiter Maximum number of CAVI iterations
#' @param elbo_tol The (percent) change in the ELBO below which it is 
#' considered converged
#' @param thin The number of iterations to wait each time before
#' re-calculating the elbo
#' @param verbose Print convergence messages
#' @param tau_q Hyperparameter tau_q
#' @param tau_mu Hyperparameter tau_mu
#' @param tau_c Hyperparameter tau_c
#' @param a Hyperparameter a
#' @param b Hyperparameter b
#' @param tau_alpha Hyperparameter tau_alpha
#' @param a_beta Hyperparameter a_beta
#' @param b_beta Hyperparameter b_beta
#' @param q Priors on the latent variables
#' @param model_mu Logical - should a gene-specific intercept term be modelled?
#' @param scale_y Logical - should the expression matrix be centre scaled?
#' @param z_init The initialisation of the latent trajectory. Should be one of
#' \enumerate{
#' \item A positive integer describing which principal component of the data should
#' be used for initialisation (default 1), \emph{or}
#' \item A numeric vector of length number of samples to be used directly for initialisation, \emph{or}
#' \item The text character \code{"random"}, for random initialisation from a standard
#' normal distribution.
#' }
#' 
#' @import Rcpp
#' @return 
#' A list whose entries correspond to the converged values of the
#' variational parameters along with the ELBO.
#' 
#' @importFrom stats coef lm prcomp
#' @importFrom stats rnorm sd
#' 
#' @useDynLib phenopath
#' 
#' @examples 
#' sim <- simulate_phenopath()
#' fit <- clvm(sim$y, sim$x)
#' 
#' 
#' @export
clvm <- function(y, x, maxiter = 1e4,
                 elbo_tol = 1e-5,
                 thin = 1,
                 verbose = TRUE,
                 z_init = 1,
                 tau_q = 1,
                 tau_mu = 1,
                 tau_c = 1,
                 a = 2, b = 2,
                 tau_alpha = 1,
                 a_beta = 10, b_beta = 1,
                 q = rep(0, nrow(y)),
                 model_mu = FALSE,
                 scale_y = TRUE) {

  N <- nrow(y)
  G <- ncol(y)
  P <- ncol(x)
  
  if(scale_y) {
    y <- scale(y)
  }
  
  ## Parameter initialisation
  
  a_tau <- rep(1, G) # rgamma(G, 2)
  b_tau <- rep(1, G) # rgamma(G, 2)


  if(P == 1) {
    f <- apply(y, 2, function(yy) {
      fit <- lm(yy ~ 0 + x[,1])
      coef(fit)
    })
    m_alpha <- matrix(f, nrow = 1, ncol = G)
  } else {
    m_alpha <- matrix(0, nrow = P, ncol = G)
  }
  s_alpha <- matrix(0.1, nrow = P, ncol = G)

  m_beta <- matrix(0, nrow = P, ncol = G)
  s_beta <- matrix(0.1, nrow = P, ncol = G)

  a_chi <- matrix(0.1, nrow = P, ncol = G)
  b_chi <- matrix(0.01, nrow = P, ncol = G)

  ## Three options for pst init:
  ## (1) A numeric vector of initial values
  ## (2) A single integer specifying which PC to initialise to
  ## (3) A text character "random" in which case they're drawn N(0,1)

  if(length(z_init) == 1 && is.numeric(z_init)) {  
    m_z <- prcomp(scale(y))$x[,z_init]
    m_z <- (m_z - mean(m_z)) / sd(m_z)
  } else if(length(z_init) == N && is.numeric(z_init)) {
    m_z <- z_init
  } else if(z_init == "random") {
    m_z <- rnorm(N)
  } else if(z_init == "pc1_with_noise") {
    m_z <- prcomp(scale(y))$x[,1]
    m_z <- (m_z - mean(m_z)) / sd(m_z)
    m_z <- rnorm(N, m_z)
  } else {
    stop("z initialisation not recognised")
  }
  
  s_z <- rep(0.1, N)

  m_lambda <- apply(y, 2, function(yy) coef(lm(yy ~ m_z))[2])
  m_lambda <- (m_lambda - mean(m_lambda)) / sd(m_lambda)
  s_lambda <- rep(0.1, G)

  m_mu <- rep(0, G)
  if(model_mu) {
    s_mu <- rep(0.1, G)
  } else {
    s_mu <- rep(0, G)
  }
  
  
  if(verbose) {
    cat("ELBO\t\tChange (%) \n")
  }

  elbo_old <- -Inf
  delta_elbo <- Inf
  elbos <- NULL
  elys <- elps <- elqs <- NULL
  ts <- cs <- betas <- taus <- chis <- alphas <- betas <- NULL
  i <- 1

  while(i < maxiter & delta_elbo > elbo_tol) {
    
    if(model_mu) {
      cumu <- cavi_update_mu(y, x, m_z, m_lambda, m_alpha, m_beta, a_tau, b_tau, tau_mu)
      m_mu <- cumu[,1]; s_mu <- cumu[,2]
    } else {
      m_mu <- rep(0,G)
      s_mu <- rep(0,G)
    }
    
    # Update lambda  
    culam <- cavi_update_lambda(y, x, m_z, s_z, m_alpha, m_beta, a_tau, b_tau,
                         m_mu, tau_c)
    m_lambda <- culam[,1]; s_lambda <- culam[,2]

    # Update tau
    cut <- cavi_update_tau(y, x, m_z, s_z, m_lambda, s_lambda, m_alpha, m_beta, s_alpha,
                           s_beta, m_mu, s_mu, a, b)
    a_tau <- cut[,1]; b_tau <- cut[,2]

    alpha_sum <- calculate_greek_sum(m_alpha, x)
    beta_sum <- calculate_greek_sum(m_beta, x)
    
    for(g in 1:G) {
      for(p in 1:P) {
        ## First calculate alpha update
        cua <- cavi_update_alpha(beta_sum, p-1, g-1, y, x, m_z, m_lambda, m_alpha, m_beta, a_tau, b_tau,
                                 m_mu, tau_alpha)
        alpha_sum <- update_greek_sum(g-1, p-1, alpha_sum, m_alpha[p,g], cua[1], x)
        m_alpha[p,g] <- cua[1] 
        s_alpha[p,g] <- cua[2]
        
        ## Calculate beta update
        cub <- cavi_update_beta(alpha_sum, p-1, g-1, y, x, m_z, s_z, m_lambda, m_alpha, m_beta, a_tau,
                                b_tau, a_chi, b_chi, m_mu)
        
        beta_sum <- update_greek_sum(g-1, p-1, beta_sum, m_beta[p,g], cub[1], x)
        m_beta[p,g] <- cub[1]
        s_beta[p,g] <- cub[2]

        # Calculate chi update
        cuch <- cavi_update_chi(m_beta[p,g], s_beta[p,g], a_beta, b_beta)
        a_chi[p,g] <- cuch[1]; b_chi[p,g] <- cuch[2] 
      }
    }

    cuz <- cavi_update_z(y, x, m_lambda, m_mu, s_lambda, m_alpha, m_beta, s_beta, a_tau, b_tau, q, tau_q)
    m_z <- cuz[,1]; s_z <- cuz[,2]

    ## calculate elbo and report
    if(i %% thin == 0) {
      elbo_vec <- calculate_elbo(y, x, m_z, s_z, m_lambda, s_lambda, m_alpha, s_alpha, m_beta, 
                             s_beta, a_tau, b_tau, a_chi, b_chi, m_mu, s_mu, q, tau_q, 
                             tau_mu, tau_c, a, b, tau_alpha, a_beta, b_beta,
                             as.integer(model_mu)) 
      
      elbo <- elbo_vec[1] + elbo_vec[2] - elbo_vec[3]

      delta_elbo <- abs((elbo - elbo_old) / elbo * 100)
      if(verbose) {
        cat(paste(elbo, delta_elbo, sep = "\t"), "\n")
      }
      elbo_old <- elbo
      elbos <- c(elbos, elbo)
    }
    i <- i + 1
  }
  
  if(i == maxiter) {
    warning("ELBO not converged")
  }

  rlist <- list(m_z = m_z, s_z = s_z, m_lambda = m_lambda, s_lambda = s_lambda, m_mu = m_mu, m_alpha = m_alpha,
                s_alpha = s_alpha, a_tau = a_tau, b_tau = b_tau,
                m_beta = m_beta, s_beta = s_beta, chi_exp = a_chi / b_chi,
                elbos = elbos)
  return(rlist)
}
  



#' Scale a vector to have mean 0 and variance 1
#' 
#' Scales vector to mean 0 variance 1 unless input standard deviation is 0
#' in which case original vector is returned
#' @keywords internal
#' @param x Input vector to scale
#' @return Scaled vector
scale_vec <- function(x) {
  std_dev <- sd(x)
  if(std_dev != 0) {
    return( (x - mean(x)) / std_dev )
  } else {
    return(x)
  }
}


