fit.pods=function(par, formula, dpar, data,  gsS, debug, hessian)
{
# Optimize log-likelihood of pod size parameters with specfied detection model
  optim(par=par, lnl.pods,formula=formula,data=data,dpar=dpar,gsS=gsS,debug=debug,hessian=hessian)
}

lnl.pods <-
function(par, formula, dpar, data,  gsS, debug)
{
# Computes negative log-likelihood for true pod size distribution
# for a given detection probability and parameters for a single observer
#
# Arguments:
#
#  par     - parameter values for pod size distribution
#  formula - the formula for the logistic detection model
#  dpar    - parameter values for detection model
#  data    - single observer data
#  gsS     - pod size calibration matrix
#  debug   - if TRUE will show iteration values
#
# Value: negative log-likelihood value at specified parameters
################################################################################
#
# Extract parameter values for detection model and set up probability matrix
# prob which has a row for each data pair and a column for each possible
# true pod size.
nmax=ncol(gsS)
fS=gammad(par,nmax)
probs=function(j)
{
  data$True=j
  dmat=model.matrix(formula,data)
  plogis(dmat%*%dpar)
}
ps=sapply(1:nmax, probs)
ps=ps/as.vector(ps%*%fS)
prob=rowSums(ps*t(fS*gsS[,data$podsize]))
neglnl=-sum(log(prob))
if(debug)
{
   cat("\npar=",par)
   cat("\nneglnl=",neglnl)
}
return(neglnl)
}

