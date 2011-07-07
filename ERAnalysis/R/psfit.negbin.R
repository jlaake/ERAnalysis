negbin.d <-
function(par,nmax,True=NULL)
{
# This function provides normalized probabilities for negative
# binomial for the range from 1 to nmax
#
# Arguments:
#
#  par  - parameters for negative binomial
#  nmax - maximum pod size
#  True - values of true pod size for a range of pod sizes (eg 4+) in which
#          a common distribution is being fitted with a relationship on
#          size with True
#
# The parameters:
#  True = NULL then log(size)=par[1], logit(p)=par[2]
#  True = NULL then log(size)=par[1]+par[2]*True, logit(p)=par[3]
  if(is.null(True))
  {
     pp=diff(c(0,pnbinom(0:(nmax-1),size=exp(par[1]),prob=plogis(par[2]))))
     ps=pp/sum(pp)
  } else
  {
     pp=do.call("rbind",lapply(True, function(tt) pnbinom(0:(nmax-1),size=exp(par[1]+par[2]*tt),prob=plogis(par[3]))))
     ps=t(apply(cbind(rep(0,nrow(pp)),pp),1,diff))/pp[,nmax]
  }
  return(ps)
}

psfit.negbin <-
function(par,psmat,True=NULL,nmax=20)
{
  ps=negbin.d(par,nmax,True)
  if(is.null(True))
    -sum(psmat%*%log(ps))
  else
    -sum(psmat*log(ps))
}

