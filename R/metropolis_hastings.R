# R/metropolis_hastings.R

#' Metropolis-Hastings Algorithm for Bayesian Inference
#'
#' This function performs the Metropolis-Hastings algorithm for Bayesian inference on the model.
#'
#' @param nsim The number of simulations to run.
#' @param burn The number of burn-in iterations.
#' @param thin The thinning interval.
#' @param L The number of levels in the model.
#' @param nL A vector of the number of samples at each level.
#'
#' @return A matrix containing the sampled values for the parameters.
#' @export
sim.theta <- function(nsim, burn, thin, L, nL) {
  set.seed(1234567)

  # Initialize matrices for parameters
  alpha = matrix(0, nrow = nsim, ncol = 3)
  alpha[1,] = c(-2, 0.5, 0.3)

  gamma = matrix(0, nrow = nsim, ncol = 3)
  gamma[1,] = c(1, 1.5, -0.2)

  PHI1 = matrix(0, nrow = nsim, ncol = 2)
  PHI1[1,] = DATA$PHI1

  PHI2 = matrix(0, nrow = nsim, ncol = 2)
  PHI2[1,] = DATA$PHI2

  beta = rep(0, nsim)
  beta[1] = beta.real

  # Matrix to store the sampled values
  theta = matrix(0, nrow = nsim, ncol = 3 + 3 + 1 + 2 + 2)
  colnames(theta) = c("beta", "alpha0", "alpha1", "alpha2", "gamma0", "gamma1", "gamma2", "PHI11", "PHI12", "PHI21", "PHI22")

  theta[1,] = c(beta[1], alpha[1,], gamma[1,], PHI1[1,], PHI2[1,])
  for(k in 2:nsim){
    #~~~~~~~~~~~~~~~~~~~~~~~~~~UPDATING beta  ~~~~~~~~~~~~~~~~~~~~~~

    yy <- abs(rnorm(1,beta[k-1],0.001))

    p1=full_beta(beta=yy,L=L0,nL=nL0,y=DATA$y,XX=DATA$XX,J=DATA$J,delta=DATA$delta,C=DATA$C,a=10,b=10,alpha=alpha[k-1,],PHI2=PHI2[k-1,])
    p2=full_beta(beta=beta[k-1],L=L0,nL=nL0,y=DATA$y,XX=DATA$XX,J=DATA$J,delta=DATA$delta,C=DATA$C,a=10,b=10,alpha=alpha[k-1,],PHI2=PHI2[k-1,])
    p2 <- ifelse(p2==0, p2+runif(1,0.00001,0.00002), p2)


    q1=pnorm(beta[k-1],yy,1)
    q2=pnorm(yy,beta[k-1],1)


    R.beta=min(1,(p1*q1)/(p2*q2))

    beta[k] <- ifelse(R.beta<runif(1), yy, beta[k-1])
    #~~~~~~~~~~~~~~~~~~~~~~~~~~UPDATING alpha  ~~~~~~~~~~~~~~~~~~~~~~

    R <- matrix(diag(c(0.01,0.01,0.01)), ncol = 3)
    yy2 <- c(rmvnorm(1, mean= (alpha[k-1,]+alpha[1,])/2, sigma=R))

    pa1=full_alpha(alpha=yy2,L=L0,nL=nL0,y=DATA$y,XX=DATA$XX,J=DATA$J,delta=DATA$delta,C=DATA$C,beta=beta[k],PHI1=PHI1[k-1,],PHI2=PHI2[k-1,],gamma=gamma[k-1,])
    pa2=full_alpha(alpha=alpha[k-1,],L=L0,nL=nL0,y=DATA$y,XX=DATA$XX,J=DATA$J,delta=DATA$delta,C=DATA$C,beta=beta[k],PHI1=PHI1[k-1,],PHI2=PHI2[k-1,],gamma=gamma[k-1,])
    pa2 <- ifelse(pa2<=0.00000001, pa2+runif(1,0.000001,0.00002), pa2)

    R.alpha=min(1,pa1/pa2)
    for(kk in 1:3){alpha[k,kk] <- ifelse(R.alpha<runif(1), yy2[kk], alpha[k-1,kk])}

    #~~~~~~~~~~~~~~~~~~~~~~~~~~UPDATING gamma  ~~~~~~~~~~~~~~~~~~~~~~


    R <- matrix(diag(c(0.01,0.01,0.01)), ncol = 3)
    yy3 <- c(rmvnorm(1, mean= (gamma[k-1,]+gamma[1,])/2, sigma=R))
    #cat("k=",k," yy3=",yy3)
    pg1=full_gamma (gamma=yy3,L=L0,nL=nL0,y=DATA$y,XX=DATA$XX,J=DATA$J,delta=DATA$delta,C=DATA$C,beta=beta[k],PHI1=PHI1[k-1,],PHI2=PHI2[k-1,],alpha=alpha[k,])
    pg2=full_gamma (gamma=gamma[k-1,],L=L0,nL=nL0,y=DATA$y,XX=DATA$XX,J=DATA$J,delta=DATA$delta,C=DATA$C,beta=beta[k],PHI1=PHI1[k-1,],PHI2=PHI2[k-1,],alpha=alpha[k,])

    pg2 <- ifelse(pg2<=0.0000001, pg2+runif(1,0.0001,0.0002), pg2)


    R.gamma=min(1,pg1/pg2)

    for(kk in 1:3){gamma[k,kk] <- ifelse(R.gamma<runif(1), yy3[kk], gamma[k-1,kk])}


    #~~~~~~~~~~~~~~~~~~~~~~~~~~UPDATING PHI1  ~~~~~~~~~~~~~~~~~~~~~~


    R <-matrix(c(0.000001,0.00000025,0.00000025,0.000001), ncol = 2)

    yy4 <- c(rmvnorm(1, mean= (PHI1[k-1,]+PHI1[1,])/2, sigma=R))


    pf1=full_PHI1(PHI1=yy4,L=L0,nL=nL0,y=DATA$y,XX=DATA$XX,J=DATA$J,delta=DATA$delta,C=DATA$C,beta=beta[k],PHI2=PHI2[k-1,],alpha=alpha[k,],
                  gamma=gamma[k,],rho=rho0,SIGMA=SIGMA0,MA=DATA$MA)
    pf2=full_PHI1(PHI1=PHI1[k-1,],L=L0,nL=nL0,y=DATA$y,XX=DATA$XX,J=DATA$J,delta=DATA$delta,C=DATA$C,beta=beta[k],PHI2=PHI2[k-1,],alpha=alpha[k,],
                  gamma=gamma[k,],rho=rho0,SIGMA=SIGMA0,MA=DATA$MA)


    pf2 <- ifelse(pf2<=0.001, pf2+runif(1,0.01,0.02), pf2)
    #pf1 <- ifelse(pf1<=0.001, pf1+runif(1,0.01,0.02), pf1)
    #cat("  pf1=",pf1,"   pf2=",pf2)
    qf1=pmvnorm(PHI1[k-1,],mean=yy4,sigma=R)
    qf2=pmvnorm(yy4,mean=PHI1[k-1,],sigma=R)
    #R.PHI1=min(1,pf1/pf2)
    R.PHI1=min(1,(pf1*qf1)/(pf2*qf2))
    #cat("  R.PHI1=",R.PHI1)


    for(kk in 1:2){PHI1[k,kk] <- ifelse(R.PHI1<runif(1), yy4[kk], PHI1[k-1,kk])}
    #~~~~~~~~~~~~~~~~~~~~~~~~~~UPDATING PHI2  ~~~~~~~~~~~~~~~~~~~~~~




    R <-matrix(c(0.000001,0.00000025,0.00000025,0.000001), ncol = 2)


    yy5 <- c(rmvnorm(1, mean=(PHI2[k-1,]+PHI2[1,])/2, sigma=R))


    pp1=full_PHI2(PHI2=yy5,L=L0,nL=nL0,y=DATA$y,XX=DATA$XX,J=DATA$J,delta=DATA$delta,C=DATA$C,beta=beta[k],PHI1=PHI1[k,],alpha=alpha[k,],
                  gamma=gamma[k,],rho=rho0,SIGMA=SIGMA0,MA=DATA$MA)

    pp2=full_PHI2(PHI2=PHI2[k-1,],L=L0,nL=nL0,y=DATA$y,XX=DATA$XX,J=DATA$J,delta=DATA$delta,C=DATA$C,beta=beta[k],PHI1=PHI1[k,],alpha=alpha[k,],
                  gamma=gamma[k-1,],rho=rho0,SIGMA=SIGMA0,MA=DATA$MA)


    pp2 <- ifelse(pp2<=0.001, pp2+runif(1,0.01,0.02), pp2)

    qp1=pmvnorm(PHI2[k-1,],mean=yy5,sigma=R)
    qp2=pmvnorm(yy5,mean=PHI2[k-1,],sigma=R)

    R.PHI2=min(1,(pp1*qp1)/(pp2*qp2))




    for(kk in 1:2){PHI2[k,kk] <- ifelse(R.PHI2<runif(1), yy5[kk], PHI2[k-1,kk])}


    #~~~~~~~~~~~~~~~~~~~~~~~~~~UPDATING tetha & result ~~~~~~~~~~~~~~~~~~~~~~

    theta[k,]=c(beta[k],alpha[k,],gamma[k,],PHI1[k,],PHI2[k,])

  }
  # Return the sampled parameters (after burn-in and thinning)
  b1 = theta[-(1:burn),]
  b2 = b1[c(seq(1, nsim-burn, thin)),]
  return(b2)
}
