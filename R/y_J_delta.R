# R/data_generation.R

#' Generate Data
#'
#' This function generates data for a Bayesian model with spatial effects and a zero-inflated discrete weibull distribution.
#'
#' @param L Number of levels.
#' @param nL A vector specifying the number of samples at each level.
#' @param beta A vector of beta values.
#'
#' @return A list containing various generated data objects including y0, y, delta, and others.
#' @export
y_J_delta <- function(L, nL, beta) {
  n = sum(nL)
  X = cbind(1, rbinom(n, 1, 0.2), rnorm(n))
  q = exp(-exp(X %*% alpha.real))
  pi = 1 / (1 + exp(-X %*% gamma.real))
  pi = sapply(1:n, function(i) max(0.01, min(0.99, pi[i])))

  y0 = rnbinom(n, mu = 4.3, size = 1.42)
  y1 = matrix(y0, nrow = 2)
  C = round(quantile(y1, probs = 0.93))

  for (jj in 1:n) {
    y1[jj] = min(C, y1[jj])
  }

  U = runif(n)
  y = matrix(ifelse(U > pi, 0, y1), nrow = L)
  delta = matrix(ifelse(y >= C, 1, 0), nrow = L)
  J = matrix(ifelse(y == 0, 1, 0), nrow = L)

  XX = array(NA, c(nL[1], ncol(X), L))
  for (ii in 1:L) {
    XX[,,ii] = X[(nL[ii]*ii-nL[ii]+1):(nL[ii]*ii),]
  }

  A = matrix(diag(2), 2, 2)
  m = apply(A, 1, sum)
  Q = diag(m) - s * A
  MA = ginv(t(Q) %*% Q) %*% t(A)
  covPHI = solve(Q) %x% SIGMA0

  PHI0 = rtmvnorm(n = 1, sigma = covPHI, lower = c(-3, -3, -3, -3), upper = c(3, 3, 3, 3))
  PHI = matrix(PHI0, ncol = 2, nrow = L, byrow = TRUE)

  PHI1 = PHI[, 1]
  PHI2 = PHI[, 2]

  return(list(n = n, m = m, MA = MA, y0 = y0, y = y, C = C, delta = delta, J = J, XX = XX,
              covPHI = covPHI, PHI0 = PHI0, PHI1 = PHI1, PHI2 = PHI2, PHI = PHI, Q = Q))
}
