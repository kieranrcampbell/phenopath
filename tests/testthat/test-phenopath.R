

context("Simulations")

test_that("simulate_phenopath() returns valid object", {
  ## Basic return structure
  expect_is(sim, "list")
  expect_equal(length(sim), 4)
  expect_identical(names(sim), c("parameters", "y", "x", "z"))
  
  ## Check everything correct size
  expect_equal(nrow(sim$y), N)
  expect_equal(ncol(sim$y), G)
  expect_equal(length(sim$x), N)
  expect_equal(length(sim$z), N)
})

test_that("simulate_phenopath() returns correct parameter specs", {
  param_df <- sim$parameters
  expect_equal(nrow(param_df), G)
  expect_identical(names(param_df), c("alpha", "lambda", "beta", "regime"))  
  
  regime_table <- table(param_df$regime)
  expect_equal(regime_table[['de']], G_de)
  expect_equal(regime_table[['pst']], G_pst)
  expect_equal(regime_table[['pst_beta']], G_pst_beta)
  expect_equal(regime_table[['de_pst_beta']], G_de_pst_beta)
})



context("phenopath() function")

test_that("phenopath() returns valid object", {
  expect_is(fit, "phenopath_fit")
  expect_equal(fit$thin, thin)
  expect_equal(fit$N, N)
  expect_equal(fit$G, G)
})
  
test_that("CAVI for CLVM has correctly sized outputs", {
  expect_equal(length(fit$m_z), N)
  expect_equal(length(fit$m_lambda), G)
  expect_equal(length(fit$m_mu), G)
  expect_equal(length(fit$a_tau), G)
  expect_equal(length(fit$b_tau), G)
  
  exp_dims <- c(1, G)
  expect_equal(dim(fit$chi_exp), exp_dims)
  expect_equal(dim(fit$m_alpha), exp_dims)
  expect_equal(dim(fit$m_beta), exp_dims)
})


test_that("phenopath() accepts ExpressionSets", {
  exprs_mat <- t(sim$y)
  pdata <- new("AnnotatedDataFrame", data.frame(x = sim$x))
  sce <- Biobase::ExpressionSet(exprs_mat, pdata)
  
  ## We'll test phenopath with three different x inputs:
  ## (1) Character vector from pData(sce)
  ## (2) Formula
  ## (3) The values themselves
  ## and check we get the same result for all
  
  set.seed(123)
  suppressWarnings(fit1 <- phenopath(sce, "x", maxiter = 4, verbose = FALSE))
  set.seed(123)
  suppressWarnings(fit2 <- phenopath(sce, ~ x, maxiter = 4, verbose = FALSE))
  set.seed(123)
  suppressWarnings(fit3 <- phenopath(sce, pData(sce)$x, maxiter = 4, verbose = FALSE))
  
  elbos <- sapply(list(fit1$elbos, fit2$elbos, fit3$elbos), tail, n = 1)
  expect_equal(length(unique(elbos)), 1)
})

