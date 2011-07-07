fit.poisson.re=function(x,re.factor,formula,sigma.formula=~1,initial,
                         maxit=1000,method="BFGS",hessian=FALSE,z=5,nmax=20)
{
# Fits a Poisson mixed-effects model to the pod size calibration data
# with random effects defined by re.factor.
#
# Arguments:
#
#  x             - dataframe
#  re.factor     - list of one or more vectors of factor variables to split x into
#                  a list of dataframes. Length of each factor variable vector must
#                  match number of rows in x.
#  formula       - formula for lambda of Poisson
#  sigma.formula - formula for sigma of normal error
#  initial       - required vector of initial parameter values
#  maxit         - maximum number of iterations (passed to optim)
#  method        - method used by optim for optimization
#  hessian       - passed to optim, if TRUE returns hessian for v-c matrix
#  z             -  defines integration bounds for normal (default is 5 for integration range of -5*sigma,5*sigma)
#  nmax          -  maximum pod size
#
#  Value: list returned with optim results
#
   mm=lapply(split(x,re.factor,drop=TRUE),function(x) model.matrix(formula,x))
   ss=lapply(split(x,re.factor,drop=TRUE),function(x) model.matrix(sigma.formula,x))
   observed.size=split(x$Estimate,re.factor,drop=TRUE)
   return(optim(initial,repois.lnl,method=method,control=list(maxit=maxit),hessian=hessian,
                 mm=mm,ss=ss,observed.size=observed.size,z=z,nmax=nmax))
}

repois.lnl=function(par,mm,ss,observed.size,z=5,nmax=20)
{
#
# Compute negative log-likelihood for a random effects Poisson model for
# the pod size calibration data. A normal 0,sigma error distribution is used
# for the random effect.  It is used for the log(lambda) where lambda is the
# Poisson intensity.
#
#  Arguments:
#   par           -  parameter values
#   mm            -  list of model matrices for lambda; list elements are split by random effect factor levels
#   ss            -  list of model matrices for sigma; list elements are split by random effect factor levels
#   observed.size -  list of vectors of observed sizes; list elements are split by random effect factor levels
#   z             -  defines integration bounds for normal (default is 5 for integration range of -5*sigma,5*sigma)
#   nmax          -  maximum pod size
#
#  Note that the mm,ss, observed.size lists need to be defined to split the data up
#  by the intersection of the factors used to define mm and ss such that for each grouping
#  there is only a single sigma.
#
#  Value:  negative log-likelihood
#
  pois.reint=function(x,obs,xbeta,sigma,nmax=20)
  {
#   Internal function to be used by integrate function for integration over normal random effect
#
#   Arguments:
#
#   x    - vector of values passed by integrate
#   obs  - observed size values
#   xbeta- x%*%beta where beta is parameter vector for lambda
#   sigma- value of sigma for normal
#   nmax - maximum pod size
#
    lambda=exp(outer(xbeta,x,"+"))
    return(dnorm(x,sd=sigma)*apply(lambda,2,function(lambda)prod(pois.d(lambda,x=obs,nmax=nmax))))
  }
# The parameter vector contains sigma parameters and then lambda parameters
# The number of columns of the design matrix ss determines number of sigma parameters
  npar=ncol(ss[[1]])
  neglnl=0
# Loop over each random effect
  for(i in 1:length(mm))
  {
#    Compute xbeta (x%*%beta) where x is design matrix and beta are parameters for lambda
     xbeta=as.vector(mm[[i]]%*%par[(npar+1):length(par)])
#    Compute sigma
     sigma=unique(exp(as.vector(ss[[i]]%*%par[1:npar])))
     if(length(sigma)>1)stop("\nInvalid sigma values\n")
     obs=observed.size[[i]]
#    Compute negative log likelihood by integrating this portion (subset of data defined by
#    random effects) of the likelihood over the normal distribution for the random effect
#    and accumulate the value.
     neglnl=neglnl-log(integrate(pois.reint,-z*sigma,z*sigma,obs=obs,xbeta=xbeta,sigma=sigma,nmax=nmax,stop.on.error=FALSE)$value)
  }
  cat("\npar = ",par)
  cat("\n neglnl= ",neglnl)
  return(neglnl)
}

fit.poisson=function(x,formula,initial,maxit=1000,method="BFGS",hessian=FALSE,nmax=20)
{
# Fits a Poisson fixed-effects model to the pod size calibration data.
#
# Arguments:
#
#  x             - dataframe
#  formula       - formula for lambda of Poisson
#  initial       - required vector of initial parameter values
#  maxit         - maximum number of iterations (passed to optim)
#  method        - method used by optim for optimization
#  hessian       - passed to optim, if TRUE returns hessian for v-c matrix
#  nmax         -  maximum pod size
#
#  Value: list returned with optim results
#
   xmat=model.matrix(formula,x)
   return(optim(initial,pois.lnl,method=method,control=list(maxit=maxit),
                hessian=hessian,xmat=xmat,observed.size=x$Estimate,nmax=nmax))
}

pois.lnl=function(par,xmat,observed.size,nmax=20)
{
# Compute negative log-likelihood for a fixed effects Poisson model for
# the pod size calibration data.
#
#  Arguments:
#   par           -  parameter values
#   xmat          -  model matrix for lambda
#   observed.size -  vector of observed sizes
#   nmax          -  maximum pod size
#
#  Value:  negative log-likelihood
#
  lambda=exp(xmat%*%par)
  neglnl=-sum(log(pois.d(lambda,x=observed.size,nmax=nmax)))
  cat("\npar = ",par)
  cat("\n neglnl= ",neglnl)
return(neglnl)
}

pois.d=function(lambda,x,nmax)
{
# Computes the Poisson probabilities for a vector of lambdas matched with observed values (x).
# Because these are pod size estimates (1,2,...) the observed x is shifted (x-1) such that it
# matches a Poisson with values 0,1,2,...  Also, to avoid dealing with an infinite number of
# of possible values the maximum pod size is set (nmax) and the distribution is renormalized
# over the range 0,1,2,...,(nmax-1).
#
# Arguments:
#   lambda        - vector of lambdas
#   x             - vector of observed pod sizes
#   nmax          -  maximum pod size
#
   if(length(x)!=length(lambda))stop("\n****invalid lengths****/n")
   num=apply(matrix(ppois(c(x-1,x-2),lambda=rep(lambda,2),lower=FALSE),ncol=2),1,diff)
   denom=ppois(rep(nmax-1,length(x)),lambda)
   denom[num<1e-16]=1
   return(num/denom)
}
