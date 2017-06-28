N <- 100
G <- 100
G_de <- 10
G_pst <- 20
G_pst_beta <- 50
G_de_pst_beta <- 20
thin <- 5


sim <- simulate_phenopath(N, G, G_de, G_pst, G_pst_beta, G_de_pst_beta)
fit <- phenopath(sim$y, sim$x, elbo_tol = 1e-3, 
                 thin = thin, verbose = FALSE)
