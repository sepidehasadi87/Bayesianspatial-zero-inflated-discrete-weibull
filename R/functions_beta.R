# R/functions_beta.R

#' Full posterior for Beta
#'
#' This function computes the full posterior probability for the Beta parameter in a Bayesian model.
#'
#' @param beta The current value of the Beta parameter.
#' @param L The number of levels in the model.
#' @param nL A vector containing the length of each level.
#' @param y The observed data.
#' @param XX The explanatory variable matrix.
#' @param J Indicator matrix for zero values in the data.
#' @param delta Indicator matrix for censored data.
#' @param C The threshold value.
#' @param a The shape parameter for the inverse gamma prior.
#' @param b The rate parameter for the inverse gamma prior.
#' @param alpha The alpha parameter in the model.
#' @param PHI2 The spatial effects parameter for PHI2.
#'
#' @return The computed posterior probability for Beta.
#' @export
full_beta <- function(beta, L, nL, y, XX, J, delta, C, a, b, alpha, PHI2) {
  wb <- matrix(NA, L, nL[1])
  for(i in 1:L) {
    for(j in 1:nL[1]) {
      w = exp(-exp((XX[j,,i]) %*% alpha + PHI2[i]))
      wb[i,j] = (((w^(y[i,j])^beta) - ((w^(y[i,j]+1))^beta)))^((1-J[i,j])*(1-delta[i,j])) *
        ((w)^C^beta)^delta[i,j]
    }
  }
  A = pinvgamma(q = beta, shape = a, rate = b, scale = 1 / b, lower.tail = TRUE, log.p = FALSE)
  p_beta = prod(wb) * A
  return(p_beta)
}
