# R/functions_alpha.R

#' Full posterior for Alpha
#'
#' This function computes the full posterior probability for the Alpha parameter in a Bayesian model.
#'
#' @param alpha The current value of the Alpha parameter.
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
#' @param gamma The Gamma parameter in the model.
#'
#' @return The computed posterior probability for Alpha.
#' @export
full_alpha <- function(alpha, L, nL, y, XX, J, delta, C, beta, PHI1, PHI2, gamma) {
  SIGMA2 = 100
  wa <- matrix(NA, L, nL[1])
  for(i in 1:L) {
    for(j in 1:nL[1]) {
      w1 <- exp(-exp((XX[j,,i]) %*% alpha + PHI2[i]))
      w2 <- exp(-XX[j,,i] %*% gamma + PHI1[i])
      wa[i,j] = (1 - (w1) + ((w2 + 1)^-1) * w1)^(J[i,j] * (1 - delta[i,j])) *
        ((w1^y[i,j]^beta) - (w1^(y[i,j] + 1)^beta))^(1 - J[i,j] * (1 - delta[i,j])) *
        (w1^C^beta)^delta[i,j]
    }
  }
  p_alpha = prod(wa) * exp(-(t(alpha) %*% ginv(SIGMA) %*% alpha) / 2)
  return(c(p_alpha))
}

