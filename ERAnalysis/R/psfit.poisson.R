psfit.poisson <-
function(par,psmat,True=NULL,nmax=20)
{
# This function computes the negative log-likelihood for fitting the Poisson
# to the calibration data
  ps=poisson.d(par,nmax,True)
  ps=log(ps)
  ps[is.infinite(ps)]=-200
  if(is.null(True))
    -sum(psmat%*%ps)
  else
    -sum(psmat*ps)
}

poisson.d <-
function(par,nmax,True=NULL)
{
# This function provides normalized probabilities for Poisson
# for the range from 1 to nmax
#
# Arguments:
#
#  par  - parameters for Poisson
#  nmax - maximum pod size
#  True - values of true pod size for a range of pod sizes (eg 4+) in which
#          a common distribution is being fitted with a relationship between
#          True and lambda
#
# The parameters:
#  True = NULL then log(lambda)=par[1]
#  True = NULL then log(lambda)=par[1]+par[2]*True

  if(is.null(True))
  {
     pp=diff(c(0,ppois(0:(nmax-1),lambda=exp(par))))
     ps=pp/sum(pp)
  }
  else
  {
     pp=do.call("rbind",lapply(True, function(tt) ppois(0:(nmax-1),lambda=exp(par[1]+par[2]*tt))))
     ps=t(apply(cbind(rep(0,nrow(pp)),pp),1,diff))/pp[,nmax]
  }
  return(ps)
}
