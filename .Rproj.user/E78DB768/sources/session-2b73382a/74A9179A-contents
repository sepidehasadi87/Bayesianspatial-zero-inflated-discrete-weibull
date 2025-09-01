# R/spatial_effects.R

#' Full posterior for PHI1
#'
#' Computes the full posterior probability for PHI1 in the model.
#'
#' @param PHI1 The current value of the PHI1 parameter.
#' @param L The number of levels.
#' @param nL The number of samples at each level.
#' @param y The observed data.
#' @param XX The explanatory variable matrix.
#' @param J Indicator matrix for zero values.
#' @param delta Indicator matrix for censored data.
#' @param C Threshold value.
#' @param beta The Beta parameter in the model.
#' @param PHI2 The PHI2 parameter.
#' @param alpha The Alpha parameter.
#' @param gamma The Gamma parameter.
#' @param rho Spatial smoothing parameter.
#' @param SIGMA The covariance matrix.
#' @param MA The matrix of neighbors.
#'
#' @return The posterior probability for PHI1.
#' @export
full_PHI1 <- function(PHI1, L, nL, y, XX, J, delta, C, beta, PHI2, alpha, gamma, rho, SIGMA, MA) {
  wp = matrix(NA, L, ncol(y))
  for (i in 1:L) {
    for (j in 1:ncol(y)) {
      w1 = exp(-exp((XX[j,,i]) %*% alpha + PHI2[i]))
      w2 = exp(-XX[j,,i] %*% gamma + PHI1[i])
      wp[i,j] = ((1 - w1 + (w2 + 1)^-1 * w1)^(J[i,j] * (1 - delta[i,j]))) *
        (1 - (w2 + 1)^-1)^(1 - J[i,j]) * (1 - (w2 + 1)^-1)^(delta[i,j] * J[i,j])
    }
  }

  B = c(PHI1 - rho * sqrt(SIGMA[1,1] / SIGMA[2,2]) * PHI2)
  D = (1 - rho^2) * SIGMA[1,1] * MA
  E = exp(-1/2 * (t(B) %*% ginv(D) %*% B))
  p_PHI1 = prod(wp) * E
  return(c(p_PHI1))
}

#' Full posterior for PHI2
#'
#' Computes the full posterior probability for PHI2 in the model.
#'
#' @param PHI2 The current value of the PHI2 parameter.
#' @param L The number of levels.
#' @param nL The number of samples at each level.
#' @param y The observed data.
#' @param XX The explanatory variable matrix.
#' @param J Indicator matrix for zero values.
#' @param delta Indicator matrix for censored data.
#' @param C Threshold value.
#' @param beta The Beta parameter in the model.
#' @param PHI1 The PHI1 parameter.
#' @param alpha The Alpha parameter.
#' @param gamma The Gamma parameter.
#' @param rho Spatial smoothing parameter.
#' @param SIGMA The covariance matrix.
#' @param MA The matrix of neighbors.
#'
#' @return The posterior probability for PHI2.
#' @export
full_PHI2 <- function(PHI2, L, nL, y, XX, J, delta, C, beta, PHI1, alpha, gamma, rho, SIGMA, MA) {
  wp = matrix(NA, L, ncol(y))
  for (i in 1:L) {
    for (j in 1:ncol(y)) {
      w1 = exp(-exp((XX[j,,i]) %*% alpha + PHI2[i]))
      w2 = exp(-XX[j,,i] %*% gamma + PHI1[i])
      wp[i,j] = ((1 - w1 + (w2 + 1)^-1 * w1)^(J[i,j] * (1 - delta[i,j]))) *
        (1 - (w2 + 1)^-1)^(1 - J[i,j]) * (1 - (w2 + 1)^-1)^(delta[i,j] * J[i,j])
    }
  }

  Sigma.p = matrix(c(1, rho, rho, 1), 2, 2)
  B = c(PHI2 - rho * sqrt(Sigma.p[2,2] / Sigma.p[1,1]) * PHI1)
  D = (1 - rho^2) * Sigma.p[2,2] * MA
  E = exp(-1/2 * (t(B) %*% ginv(D) %*% B))
  p_PHI2 = prod(wp) * E
  return(c(p_PHI2))
}

