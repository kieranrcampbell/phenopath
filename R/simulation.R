

simulate_one_gene <- function(N, pst, x, alpha = 0, c = 0, beta = 0, tau = 1e6) {
  mu <- alpha * x + (c + beta * x) * pst
  rnorm(N, mu, 1 / sqrt(tau))
}

sample_de <- function() {
  alpha <- sample(c(-1, 1), 1)
  c <- beta <- 0
  return(c(alpha, c, beta))
}

sample_pst <- function() {
  # c <- 1 * sample(c(-1, 1), 1)
  c <- rnorm(1)
  alpha <- beta <- 0
  return(c(alpha, c, beta))
}

sample_pst_beta <- function() {
  # c <- 1 * sample(c(-1, 1), 1)  
  c <- rnorm(1)
  beta <- sample(c(-1, 1), 1)  
  alpha <- 0
  return(c(alpha, c, beta))
}

sample_de_pst_beta <- function() {
  # c <- 1 * sample(c(-1, 1), 1)  
  c <- rnorm(1)
  beta <- sample(c(-1, 1), 1)  
  alpha <- sample(c(-1, 1), 1)
  return(c(alpha, c, beta))
}


#' Simulate from phenopath model
#' 
#' @export
#' @importFrom dplyr mutate
#' @importFrom tibble as_data_frame
simulate_phenopath <- function(N = 100, G = 40) {

  x <- rep(c(1, -1), each = N / 2)
  pst <- rnorm(N)
  pst <- (pst - mean(pst)) / sd(pst)

  de_pars <- t(replicate(G / 4, sample_de()))
  pst_pars <- t(replicate(G / 4, sample_pst()))
  pst_beta_pars <- t(replicate(G / 4, sample_pst_beta()))
  de_pst_beta_pars <- t(replicate(G / 4, sample_de_pst_beta()))
  
  regimes <- c("de", "pst", "pst_beta", "de_pst_beta")
  
  pars <- as_data_frame(rbind(de_pars, pst_pars, pst_beta_pars, de_pst_beta_pars))
  names(pars) <- c("alpha", "c", "beta")
  pars <- mutate(pars, regime = rep(regimes, each = G / 4))
  
  pars$alpha <- scale_vec(pars$alpha)
  pars$c <- scale_vec(pars$c)
  pars$beta <- scale_vec(pars$beta)
  
  gex <- apply(pars, 1, function(p) {
    p <- as.numeric(p)
    simulate_one_gene(N, pst, x, p[1], p[2], p[3])
  })
  
  return(list(parameters = pars, y = gex, x = x, z = pst))
}