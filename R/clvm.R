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
#' @param pst_init KIERAN WRITE
#' 
#' @return 
#' A list whose entries correspond to the converged values of the
#' variational parameters along with the ELBO.
#' 
#' @useDynLib clvm
#' 
#' 
#' @export
clvm <- function(y, x, maxiter = 1e4,
                 elbo_tol = 0.005,
                 thin = 1,
                 verbose = TRUE,
                 pst_init = 1,
                 tau_q = 1,
                 tau_mu = 1,
                 tau_c = 1,
                 a = 2, b = 2,
                 tau_alpha = 1,
                 a_beta = 10, b_beta = 1,
                 q = rep(0, nrow(y)),
                 model_mu = FALSE,
                 true_beta = NULL, true_alpha = NULL,
                 true_chi = NULL,
                 true_c = NULL,
                 true_tau = NULL,
                 scale_y = TRUE,
                 update_pst = TRUE) {

  N <- nrow(y)
  G <- ncol(y)
  P <- ncol(x)
  
  if(scale_y) {
    y <- scale(y, scale = T)
  }
  
  ## Parameter initialisation
  
  if(is.null(true_tau)) {
    a_tau <- rep(1, G) # rgamma(G, 2)
    b_tau <- rep(1, G) # rgamma(G, 2)
  } else {
    a_tau = true_tau^2 / 0.001
    b_tau = true_tau / 0.001
  }

  if(is.null(true_alpha)) {
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
  } else {
    m_alpha = true_alpha
    s_alpha <- matrix(0.001, nrow = P, ncol = G)
  }
  
  
  if(is.null(true_beta)) {
    m_beta <- matrix(0, nrow = P, ncol = G)
    s_beta <- matrix(0.1, nrow = P, ncol = G)
  } else {
    m_beta = true_beta
    s_beta <- matrix(0.0000001, nrow = P, ncol = G)
  }
  
  if(is.null(true_chi)) {
    a_chi <- matrix(0.1, nrow = P, ncol = G)
    b_chi <- matrix(0.01, nrow = P, ncol = G)
  } else {
    a_chi <- matrix(true_chi[,1], nrow = 1)
    b_chi <- matrix(true_chi[,2], nrow = 1)
  }
  
  ## Three options for pst init:
  ## (1) A numeric vector of initial values
  ## (2) A single integer specifying which PC to initialise to
  ## (3) A text character "random" in which case they're drawn N(0,1)

  if(length(pst_init) == 1 && is.numeric(pst_init)) {  
    m_t <- prcomp(scale(y))$x[,pst_init]
    m_t <- (m_t - mean(m_t)) / sd(m_t)
  } else if(length(pst_init) == N && is.numeric(pst_init)) {
    m_t <- pst_init
  } else if(pst_init == "random") {
    m_t <- rnorm(N)
  } else if(pst_init == "pc1_with_noise") {
    m_t <- prcomp(scale(y))$x[,1]
    m_t <- (m_t - mean(m_t)) / sd(m_t)
    m_t <- rnorm(N, m_t)
  } else {
    stop("z initialisation not recognised")
  }
  
  if(update_pst == TRUE) {
    s_t <- rep(0.1, N)
  } else {
    s_t <- rep(0.00001, N)
  }
  
  
  if(is.null(true_c)) {
    m_c <- apply(y, 2, function(yy) coef(lm(yy ~ m_t))[2])
    m_c <- (m_c - mean(m_c)) / sd(m_c)
    s_c <- rep(0.1, G)
  } else {
    m_c = true_c
    s_c <- rep(0.0001, G)
  }
  
  m_mu <- rep(0, G)
  s_mu <- rep(0.1, G)
  
  
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
      cumu <- cavi_update_mu(y, x, m_t, m_c, m_alpha, m_beta, a_tau, b_tau, tau_mu)
      m_mu <- cumu[,1]; s_mu <- cumu[,2]
    } else {
      m_mu <- rep(0,G)
      s_mu <- rep(0,G)
    }
    
    cuc <- cavi_update_c(y, x, m_t, s_t, m_alpha, m_beta, a_tau, b_tau,
                         m_mu, tau_c)
    if(is.null(true_c)) {
      m_c <- cuc[,1]; s_c <- cuc[,2]
    }
    cs <- c(cs, sum(m_c^2))
    
    cut <- cavi_update_tau(y, x, m_t, s_t, m_c, s_c, m_alpha, m_beta, s_alpha,
                           s_beta, m_mu, s_mu, a, b)
    
    if(is.null(true_tau)) {
      a_tau <- cut[,1]; b_tau <- cut[,2]
    }
    taus <- c(taus, sum((a_tau / b_tau)^2))
    
    alpha_sum <- calculate_greek_sum(m_alpha, x)
    beta_sum <- calculate_greek_sum(m_beta, x)
    
    for(g in sample(1:G)) {
      for(p in sample(1:P)) {
        ## First calculate alpha update
        cua <- cavi_update_alpha(beta_sum, p-1, g-1, y, x, m_t, m_c, m_alpha, m_beta, a_tau, b_tau,
                                 m_mu, tau_alpha)
        
        ## Now sort the modified _alpha_ sum matrix
        # for(i in 1:N) { # bite me
        #   alpha_sum[g,i] <- -m_alpha[p,g] * x[i,p] + cua[1] * x[i,p]
        # }
        
        
        
        # Now update the alphas
        if(is.null(true_alpha)) {
          alpha_sum <- update_greek_sum(g-1, p-1, alpha_sum, m_alpha, cua[1], x)
          m_alpha[p,g] <- cua[1] 
          s_alpha[p,g] <- cua[2]
        }
        
        ## Calculate beta update
        cub <- cavi_update_beta(alpha_sum, p-1, g-1, y, x, m_t, s_t, m_c, m_alpha, m_beta, a_tau,
                                b_tau, a_chi, b_chi, m_mu)
        
        # now sort the betas
        # for(i in 1:N) {
        #   beta_sum[g,i] <- -m_beta[p,g] * x[i,p] + cub[1] * x[i,p]
        # }
        
        if(is.null(true_beta)) {
          beta_sum <- update_greek_sum(g-1, p-1, beta_sum, m_beta, cub[1], x)
          m_beta[p,g] <- cub[1]
          s_beta[p,g] <- cub[2]
        }
                
        cuch <- cavi_update_chi(m_beta[p,g], s_beta[p,g], a_beta, b_beta)
        if(is.null(true_chi)) {
          a_chi[p,g] <- cuch[1]; b_chi[p,g] <- cuch[2] 
        }
      }
    }
    
    alphas <- c(alphas, sum(m_alpha^2))
    betas <- c(betas, sum(m_beta^2))

    cup <- cavi_update_pst(y, x, m_c, m_mu, s_c, m_alpha, m_beta, s_beta, a_tau, b_tau, q, tau_q)
    if(update_pst) {
      m_t <- cup[,1]; s_t <- cup[,2]
    }
    ts <- c(ts, sum(m_t^2))
    
    ## calculate elbo and report
    if(i %% thin == 0) {
      elbo_vec <- calculate_elbo(y, x, m_t, s_t, m_c, s_c, m_alpha, s_alpha, m_beta, 
                             s_beta, a_tau, b_tau, a_chi, b_chi, m_mu, s_mu, q, tau_q, 
                             tau_mu, tau_c, a, b, tau_alpha, a_beta, b_beta,
                             as.integer(model_mu)) 
      
      # print(elbo_vec)
      elbo <- elbo_vec[1] + elbo_vec[2] - elbo_vec[3]
      elys <- c(elys, elbo_vec[1])
      elps <- c(elps, elbo_vec[2])
      elqs <- c(elqs, elbo_vec[3])
      
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

  rlist <- list(m_t = m_t, s_t = s_t, m_c = m_c, s_c = s_c, m_mu = m_mu, m_alpha = m_alpha,
                s_alpha = s_alpha, a_tau = a_tau, b_tau = b_tau,
                m_beta = m_beta, s_beta = s_beta, chi_exp = a_chi / b_chi,
                elbos = elbos, elys = elys, elps = elps, elqs = elqs,
                cs = cs, taus = taus, ts = ts, alphas = alphas, betas = betas)
  return(rlist)
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
significant_interactions <- function(pcavi, n = 2) {
  m_beta <- pcavi$m_beta
  pos_sd <- sqrt(pcavi$s_beta)
  
  sig <- sapply(seq_len(nrow(m_beta)), function(i) {
    m_beta[i,] - n * pos_sd[i,] > 0 | m_beta[i,] + n * pos_sd[i,] < 0
  })
  return(sig)
}

#' Plots the ELBO
#' @export
#' @param fit A CLVM fit
#' @import ggplot2
plot_elbo <- function(fit) {
  elbo <- fit$elbos[-(1)]
  qplot(seq_along(elbo) + 1, elbo, geom = c("point", "line")) + xlab("Iter") + ylab("ELBO")
}


