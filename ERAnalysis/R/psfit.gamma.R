psfit.gamma <-  function(par,psmat,True=NULL,shape=TRUE,nmax=20)
{
# This function computes the negative log-likelihood for fitting the negative
# binomial to the calibration data  ps=gammad(par,20,True,shape)
  ps=gammad(par,nmax,True,shape)
  if(is.null(True))
    -sum(psmat%*%log(ps))
  else
    -sum(psmat*log(ps))
}

gammad <- function(par,nmax,True=NULL,shape=TRUE)
{
# This function creates a discretized gamma function from 1 to nmax
#
# Arguments:
#
#  par  - parameters for gamma
#  nmax - maximum pod size
#  True - values of true pod size for a range of pod sizes (eg 4+) in which
#          a common distribution is being fitted with a relationship on
#          either the shape (if shape=TRUE), or the scale=1/rate
#  shape - TRUE/FALSE; see above
#
# The parameters:
#  True = NULL then log(shape)=par[1], log(rate)=par[2] and scale=1/rate
#  True is not null and !shape then log(shape)=par[1], log(rate)=-par[2]-par[3]*True
#  True is not null and shape=TRUE then log(shape)=par[2]+par[3]*True, log(rate)=par[1]
  if(is.null(True))
  {
     pp=pgamma(1:nmax,shape=exp(par[1]),rate=exp(par[2]))
     ps=diff(c(0,pp))/pp[nmax]
     ps[ps==0]=1e-12
     ps=ps/sum(ps)
  }
  else
  {
     if(!shape)
     {
        pp= do.call("rbind",lapply(True, function(x) pgamma(1:nmax,shape=exp(par[1]),rate=exp(-par[2]-par[3]*x))))
        ps=t(apply(cbind(rep(0,nrow(pp)),pp),1,diff))/pp[,nmax]
     }
     else
     {
        pp= do.call("rbind",lapply(True, function(x) pgamma(1:nmax,shape=exp(par[2]+par[3]*x),rate=exp(par[1]))))
        ps=t(apply(cbind(rep(0,nrow(pp)),pp),1,diff))/pp[,nmax]
     }
     ps[ps==0]=1e-12
     ps=ps/rowSums(ps)
  }
  return(ps)
}
