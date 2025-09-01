# R/functions_gamma.R

#' Full posterior for Gamma
#'
#' This function computes the full posterior probability for the Gamma parameter in a Bayesian model.
#'
#' @param gamma The current value of the Gamma parameter.
#' @param L The number of levels in the model.
#' @param nL A vector containing the length of each level.
#' @param y The observed data.
#' @param XX The explanatory variable matrix.
#' @param J Indicator matrix for zero values in the data.
#' @param delta Indicator matrix for censored data.
#' @param C The threshold value.
#' @param beta The Beta parameter in the model.
#' @param PHI1 The spatial effects parameter for PHI1.
#' @param PHI2 The spatial effects parameter for PHI2.
#' @param alpha The Alpha parameter in the model.
#'
#' @return The computed posterior probability for Gamma.
#' @export
full_gamma <- function(gamma, L, nL, y, XX, J, delta, C, beta, PHI1, PHI2, alpha) {
  SIGMA2 = 100
  wg <- matrix(NA, L, nL[1])
  for(i in 1:L) {
    for(j in 1:nL[1]) {
      w1 <- exp(-exp((-XX[j,,i] %*% alpha + PHI2[i])))
      w2 <- exp(-XX[j,,i] %*% gamma + PHI1[i])
      wg[i,j] = (1 - w1 + (w2 + 1)^-1 * w1)^(J[i,j] * (1 - delta[i,j])) *
        (1 - (w2 + 1)^-1)^(1 - J[i,j]) * (1 - (w2 + 1)^-1)^(delta[i,j] * J[i,j])
    }
  }
  p_gamma = prod(wg) * exp((-t(gamma) %*% ginv(SIGMA) %*% gamma) / 2)
  return(c(p_gamma))
}
