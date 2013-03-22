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
###  Inputs:      pred = a vector or design matrix of predictor values (i.e., counterfactuals)
##           n.sims.fe = number of fixed effect draws  (default = 10,000)
##           n.sims.re = number of random effect draws (maximum/default = 10,000)
##             coef.fe = mean estimates of fixed effects
##             vcov.fe = the variance-covariance matrix of the fixed effects
##             vcov.re = the variance-covariance matrix of the random effects
##                link = link function ("logit", "log", default: "identity")
##                  ci = confidence interval (default = 0.95)
##
##   To-do:  - Extend to 3-level or higher designs
##           - Computing first differences
##
### Reference:
##    King, G., Tomz, M., & Wittenberg, J. (2000). Making the most of statistical
##    analyses: Improving interpretation and presentation.
##    American Journal of Political Science, 44, 347–361. doi:10.3886/ICPSR01255.v1
##

sim.marginal.pred <- function(pred, n.sims.fe=10000, n.sims.re=10000, coef.fe, vcov.fe, vcov.re, link="identity", ci=0.95) {
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
  
  # simulate fixed effect draws
  sim.fe <- mvrnorm(n=n.sims.fe, mu=coef.fe, Sigma=vcov.fe)

  # cap the max number of random effects at 10,000 (so R isn't overwhelmed)
  num.re <- ncol(vcov.re)  # number of random effects
  n.sims.re <- min(c(n.sims.re, 10000))
  
  ## Generate predictions across covariate combinations
  pred.mat <- t(rbind(pred))   # coerce covariate values to rows
  num.pred <- ncol(pred.mat)   # determine number of covariate combos

  # empty matrix to collect simulated results
  simyn <- matrix(NA, nrow=n.sims.fe, ncol=num.pred)
  
  for (j in 1:n.sims.fe) {
    sim.re <- mvrnorm(n=n.sims.re, mu=rep(0, num.re), Sigma=vcov.re)
  
    simmu.fe <- sim.fe[j,] %*% pred.mat
    simmu.re <- sim.re %*% rep(1, num.re)
    simmu <- apply(simmu.fe, 2, function(x, re) x + re, re=simmu.re[,1])
    simy0 <- unlink(simmu, link)      # create vector of simulated predictions with diff't random effects
    simyn[j,] <- apply(simy0, 2, mean) # average over random effects
  }
    
  # output mean, and 95% CI limits
  simci <- HPDinterval(as.mcmc(simyn), prob=ci)
  simmean <- apply(simyn, 2, mean)
  pred.out <- cbind(simmean, simci[,"lower"], simci[,"upper"])
  colnames(pred.out) <- c("mean","lower","upper")
  
  # Return predicted values
  return(pred.out)
}

