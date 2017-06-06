
#' Plots the ELBO
#' @export
#' @param fit A CLVM fit
#' @import ggplot2
plot_elbo <- function(fit) {
  elbo <- fit$elbos[-(1)]
  qplot(seq_along(elbo) + 1, elbo, geom = c("point", "line")) + xlab("Iter") + ylab("ELBO")
}
