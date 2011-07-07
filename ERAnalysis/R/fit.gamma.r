fit.gamma.re=function(x,re.factor,sformula,rformula,sigma.formula=~1,initial,maxit=1000,method="BFGS",hessian=FALSE,nmax=20,z=5)
{
# Fits a gamma mixed-effects model to the pod size calibration data.
#
# Arguments:
#
#  x             - dataframe
#  re.factor     - list of one or more vectors of factor variables to split x into
#                  a list of dataframes. Length of each factor variable vector must
#                  match number of rows in x.
#  sformula      - formula for shape of gamma
#  rformula      - formula for rate of gamma
#  sigma.formula - formula for sigma of normal error
#  initial       - required vector of initial parameter values
#  maxit         - maximum number of iterations (passed to optim)
#  method        - method used by optim for optimization
#  hessian       - passed to optim, if TRUE returns hessian for v-c matrix
#  nmax         -  maximum pod size
#
#  Value: list returned with optim results
#
   mm=lapply(split(x,re.factor,drop=TRUE),function(x) model.matrix(rformula,x))
   rr=lapply(split(x,re.factor,drop=TRUE),function(x) model.matrix(sformula,x))
   ss=lapply(split(x,re.factor,drop=TRUE),function(x) model.matrix(sigma.formula,x))
   observed.size=split(x$Estimate,re.factor,drop=TRUE)
   return(optim(initial,regam.lnl,method=method,control=list(maxit=maxit),
                hessian=hessian,mm=mm,rr=rr,ss=ss,observed.size=observed.size,nmax=nmax,z=z))
}

regam.lnl=function(par,mm,rr,ss,observed.size,nmax=20,z=5)
{
#
# Compute negative log-likelihood for a random effects gamma model for
# the pod size calibration data. A normal 0,sigma error distribution is used
# for the random effect.  It is used for the log(rate) where rate is a
# gamma parameter.
#
#  Arguments:
#   par           -  parameter values
#   mm            -  list of model matrices for rate; list elements are split by random effect factor levels
#   rr            -  list of model matrices for shape; list elements are split by random effect factor levels
#   ss            -  list of model matrices for sigma; list elements are split by random effect factor levels
#   observed.size -  list of vectors of observed sizes; list elements are split by random effect factor levels
#   z             -  defines integration bounds for normal (default is 5 for integration range of -5*sigma,5*sigma)
#   nmax          -  maximum pod size
#
#  Note that the mm,rr,ss, observed.size lists need to be defined to split the data up
#  by the intersection of the factors used to define mm,rr and ss such that for each grouping
#  there is only a single sigma.
#
#  Value:  negative log-likelihood
#
  gam.reint=function(x,obs,xbeta,shape,sigma,nmax=20)
  {
#   Internal function to be used by integrate function for integration over normal random effect
#
#   Arguments:
#
#   x    - vector of values passed by integrate
#   obs  - observed size values
#   xbeta- x%*%beta where beta is parameter vector for rate of gamma
#   shape- shape parameters for gamma
#   sigma- value of sigma for normal
#   nmax - maximum pod size
#
    rate=exp(outer(xbeta,x,"+"))
    return(dnorm(x,sd=sigma)*apply(rate,2,function(rate) prod(gam.d(shape,rate,x=obs,nmax=nmax))))
  }
# The parameter vector contains sigma parameters, rate parameters and then shape parameters
# The number of columns of the design matrix ss determines number of sigma parameters
  npars=ncol(ss[[1]])
  nparr=ncol(mm[[1]])
  neglnl=0
  cat("\npar = ",par)
# Loop over each random effect
  for(i in 1:length(mm))
  {
#    Compute xbeta (x%*%beta) where x is design matrix and beta are parameters for rate
     xbeta=as.vector(mm[[i]]%*%par[(npars+1):(npars+nparr)])
#    Compute sigma
     sigma=unique(exp(as.vector(ss[[i]]%*%par[1:npars])))
#    Compute shape
     shape=unique(exp(as.vector(rr[[i]]%*%par[(npars+nparr+1):length(par)])))
     if(length(sigma)>1)stop("\nInvalid sigma values\n")
     obs=observed.size[[i]]
#    Compute negative log likelihood by integrating this portion (subset of data defined by
#    random effects) of the likelihood over the normal distribution for the random effect
#    and accumulate the value.
     neglnl=neglnl-log(integrate(gam.reint,-z*sigma,z*sigma,obs=obs,xbeta=xbeta,sigma=sigma,shape=shape,nmax=nmax,stop.on.error=FALSE)$value)
  }
  cat("\npar = ",par)
  cat("\n neglnl= ",neglnl)
  return(neglnl)
}

fit.gamma=function(x,sformula,rformula,initial,maxit=1000,method="BFGS",hessian=FALSE,nmax=20)
{
# Fits a gamma fixed-effects model to the pod size calibration data.
#
# Arguments:
#
#  x             - dataframe
#  sformula      - formula for shape of gamma
#  rformula      - formula for rate of gamma
#  initial       - required vector of initial parameter values
#  maxit         - maximum number of iterations (passed to optim)
#  method        - method used by optim for optimization
#  hessian       - passed to optim, if TRUE returns hessian for v-c matrix
#  nmax         -  maximum pod size
#
#  Value: list returned with optim results
#
   smat=model.matrix(sformula,x)
   rmat=model.matrix(rformula,x)
   return(optim(initial,gam.lnl,method=method,control=list(maxit=maxit),
                hessian=hessian,rmat=rmat,smat=smat,observed.size=x$Estimate,nmax=nmax))
}

gam.lnl=function(par,smat,rmat,observed.size,nmax=20)
{
# Compute negative log-likelihood for a fixed effects gamma model for
# the pod size calibration data.
#
#  Arguments:
#   par           -  parameter values
#   smat          -  model matrix for shape
#   rmat          -  model matrix for rate
#   observed.size -  vector of observed sizes
#   nmax          -  maximum pod size
#
#  Value:  negative log-likelihood
#
  shape=exp(smat%*%par[1:ncol(smat)])
  rate=exp(rmat%*%par[(ncol(smat)+1):length(par)])
  neglnl=-sum(log(gam.d(shape,rate,x=observed.size,nmax=nmax)))
  cat("\npar = ",par)
  cat("\n neglnl= ",neglnl)
return(neglnl)
}

gam.d=function(shape,rate,x,nmax)
{
# Computes the gamma probabilities for a vector of shape and rates matched to observed values (x).
# Because these are pod size estimates (1,2,...) the integral for x is from x-1 to x of the gamma.
# Also, to avoid dealing with an infinite number of possible values the maximum pod size
# is set (nmax) and the distribution is renormalized over the range 0 to nmax.
#
# Arguments:
#   shape        - vector of shapes
#   rate         - vector of rates
#   x            - vector of observed pod sizes
#   nmax         -  maximum pod size
#
   if(length(x)!=length(rate))stop("\n****invalid lengths****/n")
   num=apply(matrix(pgamma(c(x-1,x),shape=shape,rate=rep(rate,2)),ncol=2),1,diff)
   denom=pgamma(rep(nmax,length(x)),shape=shape,rate=rate)
#   denom[num<1e-16]=1
   return(num/denom)
}

