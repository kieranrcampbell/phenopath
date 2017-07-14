
#' Plots the ELBO
#' 
#' Plots the evidence lower bound (ELBO) as a function of iterations
#' @export
#' @param fit An object returned by a call to \code{phenopath}
#' @import ggplot2
#' @examples 
#' sim <- simulate_phenopath()
#' fit <- phenopath(sim$y, sim$x)
#' plot_elbo(fit)
#' @return A \code{ggplot2} object of the ELBO against the number of iterations
plot_elbo <- function(fit) {
  stopifnot(is(fit, "phenopath_fit"))
  elbo <- fit$elbos[-(1)]
  ggplot2::qplot(fit$thin * seq_along(elbo) + 1, elbo, geom = c("point", "line")) + xlab("Iter") + ylab("ELBO")
}
