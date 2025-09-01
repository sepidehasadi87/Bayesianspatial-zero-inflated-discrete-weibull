#' Full posterior for SIGMA
#'
#' Computes the full posterior for the SIGMA parameter in the model.
#'
#' @param PHI The spatial effects parameter for PHI.
#' @param L The number of levels in the model.
#' @param rho Spatial smoothing parameter.
#' @param nu0 Degrees of freedom for the Wishart prior.
#' @param S0 Scale matrix for the Wishart prior.
#' @param sigma The covariance matrix.
#'
#' @return The computed posterior for SIGMA.
#' @export
full_SIGMA <- function(PHI, L, rho, nu0, S0, sigma) {
  SQ = t(PHI) %*% sigma %*% PHI
  SIGMA.p = riwish((nu0 + L - 1), S0 + SQ)
  return(SIGMA.p)
}


