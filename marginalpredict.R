### Description:
##    Simulates predicted values with confidence intervals for fixed and
##    random effect regression models using the multivariate normal approach
##    suggested by King, Tomz, and Wittenberg (2000).
##
##    These are marginal predictions that average over the random effects.
##   
##    Random effects models with up to 3-levels of nesting are permitted.
##    Since the algorithm performs a nested averaging, adding a third level
##    increases the computational demands exponentially (i.e., days).
##    First run the simulation with fewer iterations to gauge estimation time.
##
##    This function currently supports models using the identity, logit, and
##    log links. This can be easily extended by updating the "unlink"
##    utility function.
##
###  Original Author: David Huh
##
###  Dependencies: coda, MASS
##
###  Arguments:  pred.fe  = a vector/design matrix of hypothetical predictor values
##               pred.re  = a vector/design matrix of the level 2 random effects values
##               pred.re2 = a vector/design matrix of the level 3 random effects values
##               coef.fe  = mean estimates of fixed effects
##               vcov.fe  = the variance-covariance matrix of the fixed effects
##               vcov.re  = a variance-covariance matrix of the level 2 random effects
##               vcov.re2 = a variance-covariance matrix of the level 3 random effects
##             n.sims.fe  = number of fixed effect draws
##             n.sims.re  = number of random effect draws of the level 2 random effects
##             n.sims.re2 = number of random effect draws of the level 3 random effects
##                  link  = link function ("logit", "log", "identity" [default])
##                    ci  = confidence interval (default = 0.95)
##
###  Values:       sims = a matrix of simulations (rows: simulates, columns: predictor combo)
##              summary = a matrix summary the mean and confidence intervals across predictors
##
##   To-do:  - Computing first differences
##           - Improve efficiency via vectorization and other optimizations
##
### Reference:
##    King, G., Tomz, M., & Wittenberg, J. (2000). Making the most of statistical
##    analyses: Improving interpretation and presentation.
##    American Journal of Political Science, 44, 347-361. doi:10.3886/ICPSR01255.v1
##

marginalpredict <- function(pred.fe, pred.re, pred.re2, coef.fe, vcov.fe, vcov.re, vcov.re2,
                            n.sims.fe=100, n.sims.re=100, n.sims.re2=100,
                            link="identity", ci=0.95) {
  require(coda, quietly=TRUE)
  require(MASS, quietly=TRUE)
  
  ## validate required arguments
  if (missing(pred))
    stop("Missing a vector or design matrix of predictor values.") 
  if (missing(coef.fe))
    stop("Missing mean estimates of the fixed effect.") 
  if (missing(vcov.fe))
    stop("Missing a variance-covariance matrix for the fixed effects.")
  if (missing(vcov.re) & !missing(vcov.re2))
    stop("Missing a variance-covariance matrix for the level 2 random effects.")
  
  ## utility to transform linear predictor to units of the outcome
  unlink <- function(mu, link) {
    if (link=="logit") {
      mu.t <- 1/(1+exp(-mu))
    } else if (link=="log") {
      mu.t <- exp(mu)
    } else mu.t <- mu
    return(mu.t)
  }
  
  ## simulate draws from the fixed effect distribution
  sim.fe <- mvrnorm(n=n.sims.fe, mu=coef.fe, Sigma=vcov.fe)
  
  ## pre-process design matrix (or vector) of predictors
  pred.fe.mat <- t(rbind(pred.fe))   # coerce covariate values to rows
  num.fe.pred <- ncol(pred.fe.mat)   # determine number of covariate combos
  
  ## initialize empty matrix to collect simulated results
  simyn <- matrix(NA, nrow=n.sims.fe, ncol=num.pred)
  
  if (missing(vcov.re)) {   ## Fixed effect model ##
    simmu <- sim.fe %*% pred.fe.mat  # multiply out FE section of linear predictor
    simyn <- unlink(simmu, link)  # convert linear predictor to units of the outcome
  } else {                  ## Random effect(s) model ##
    ## pre-process design matrix (or vector) of predictors
    pred.re.mat <- t(rbind(pred.re))   # coerce covariate values to rows
    
    ## initialize variables for collecting random effects simulations
    n.sims.re <- min(c(n.sims.re, 10000)) # limit the max number of random effects
    num.re <- ncol(vcov.re)               # number of level 2 random effects
    
    if (!missing(vcov.re2)) {
      n.sims.re2 <- min(c(n.sims.re2, 10000))                 # limit the max number of random effects
      num.re2 <- ncol(vcov.re2)                               # number of level 3 random effects
      simynr <- matrix(NA, nrow=n.sims.re, ncol=num.fe.pred)  # initialize empty matrix to collect results
    }
    
    for (j in 1:n.sims.fe) {
      simmu.fe <- sim.fe[j,] %*% pred.fe.mat    # multiply out FE section of linear predictor
      
      sim.re <- mvrnorm(n=n.sims.re, mu=rep(0, num.re), Sigma=vcov.re)
      simmu.re <- sim.re %*% pred.re.mat  # multiply out RE section of linear predictor
      simmu <- matrix(rep(simmu.fe, each=n.sims.re), nrow=n.sims.re, byrow=TRUE) + simmu.re  # assemble linear predictor
      
      if (missing(vcov.re2)) {
        simy0 <- unlink(simmu, link)    # convert linear predictor to units of the outcome
        simyn[j,] <- colMeans(simy0)    # average over level 2 random effects
      } else {  ## Average over add'l level of random effects ##
        for (k in 1:n.sims.re) {
          sim.re2 <- mvrnorm(n=n.sims.re2, mu=rep(0, num.re2), Sigma=vcov.re2)
          simmu.re2 <- sim.re2 %*% pred.re2.mat  # multiply out RE section of linear predictor
          simmu2 <- matrix(rep(simmu[k,], each=n.sims.re2), nrow=n.sims.re2, byrow=TRUE) + simmu.re2  # assemble linear predictor
          
          simy0 <- unlink(simmu2, link)    # convert linear predictor to units of the outcome
          simynr[k,] <- colMeans(simy0)    # average over level 3 random effects
        }
        simyn[j,] <- colMeans(simynr)      # average over level 2 random effects
      }
    }
  }
  
  ## calculate mean and confidence intervals
  simci <- HPDinterval(as.mcmc(simyn), prob=ci)
  simmean <- colMeans(simyn)
  pred.out <- cbind(simmean, simci[,"lower"], simci[,"upper"])
  colnames(pred.out) <- c("mean","lower","upper")
  rownames(pred.out) <- seq(1, num.pred)
  
  ## return predicted values and marginal simulates
  res <- list(sims=simyn, summary=pred.out)
  return(res)
}
