### Description:
##    This function simulates predicted values w/ confidence intervals for 
##    mixed effect regression models using the multivariate normal approach
##    suggested by King, Tomz, and Wittenberg (2000).
##
##    These are marginal predictions that average over the random effects.
##   
##    This current version works for 2-level models.  For additional levels
##    of nesting, it may be necessary to create a nested averaging of random
##    effects at each level.
##
##    This model currently supports models using the identity, logit, and
##    log links. This can be easily extended by updating the "unlink"
##    utility function.
##
###  Original Author: David Huh
##
###  Dependencies: MASS, coda
##
###  Inputs:     sim.fe = a matrix of simulated fixed effect coefficients
##             pred.vec = a vector of predictor values (i.e., counterfactuals)
##                    S = the variance-covariance matrix of the random effects
##
### Reference:
##    King, G., Tomz, M., & Wittenberg, J. (2000). Making the most of statistical
##    analyses: Improving interpretation and presentation.
##    American Journal of Political Science, 44, 347–361. doi:10.3886/ICPSR01255.v1
##
##

sim.marginal.pred <- function(sim.fe, pred, S, link="identity", ci=0.95) {
  require(coda, quietly=TRUE)
  require(MASS, quietly=TRUE)
  
  unlink <- function(mu,link) {
    if (link=="logit") {
      mu.t <- 1/(1+exp(-mu))
    } else if (link=="log") {
      mu.t <- exp(mu)
    } else mu.t <- mu
    return(mu.t)
  }
  
  # simulate same number of random effects as fixed effects, but cap the max
  # number of random effects 10,000 (so R isn't overwhelmed)
  n.sims.fe <- nrow(sim.fe)
  n.sims.re <- min(c(10000, n.sims.fe))
  
  num.re <- ncol(S)  # number of random effects
  sim.re <- mvrnorm(n=n.sims.re, mu=rep(0, num.re), Sigma=S)
  simmu.re <- sim.re %*% rep(1, num.re)
  
  # empty vector to collect simulated results
  simyn <- rep(NA, n.sims.fe)
  
  for (i in 1:n.sims.fe) {
    simmu.fe <- sim.fe[i,] %*% pred
    simmu <- simmu.fe + simmu.re[,1]
    simy0 <- unlink(simmu, link)  # create vector of simulated predictions with diff't random effects
    simyn[i] <- mean(simy0)       # average over random effects
  }
  
  simci <- HPDinterval(as.mcmc(simyn), prob=ci)
  
  # output mean, and 95% CI limits
  res <- c(mean(simyn), simci[,"lower"], simci[,"upper"])
  return(res)
}
