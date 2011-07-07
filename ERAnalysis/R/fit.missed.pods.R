lnl.missed.pods <-
function(par,data,primary.only,years=NULL,formula,gsS,debug)
{
# Computes negative log-likelihood for true pod size distribution and
# detection probability parameters for double observer survey data
#
# Arguments:
#
#  par     - parameter values
#  data    - match data ordered in pairs with station=P then S where
#            P=Primary & S=Secondary
#  primary.only - data from observations of primary observer only on watch
#  years   - vector of years included in the analysis
#  formula - the formula for the logistic detection model
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
pbyyear=FALSE
if(is.null(years)) 
  npar=2
else
{
  pbyyear=TRUE
  npar=2*length(levels(factor(years)))
}
beta=par[(npar+1):length(par)]
n=nrow(data)
# Loop over each possible true pod size and compute the probability of observing
# the value of the capture history if the true size was i for each given
# station's recorded pod size if seen by the observer at the station.
# See pdf writeup.
prob=sapply(1:nmax,function(x) 
{
  data$True=x
  dmat=model.matrix(formula,data)
  p=plogis(dmat%*%beta)
  gs=(gsS[x,][data$podsize])^data$seen
  pp=p^data$seen*(1-p)^(1-data$seen)
  return( gs[seq(1,n,2)]*gs[seq(2,n,2)]*pp[seq(1,n,2)]*pp[seq(2,n,2)]/
               (1-(1-p[seq(1,n,2)])*(1-p[seq(2,n,2)])))
})
# Next compute the overall probability by summing across all true pod sizes
# weighted by the estimated probability that a pod is of that true size (fS).
# The negative sum of the log(total prob) is the negative log-likelihood.
# If pbyyear is TRUE then this is done by survey year because the fS values
# differ by year.
if(!pbyyear)
{
  fS=gammad(par[1:2],nmax)
  neglnl=-sum(log(prob%*%fS))
  if(!is.null(primary.only))neglnl=neglnl+ lnl.pods(par[1:2], formula=formula, dpar=beta, 
                              data=primary.only, gsS=gsS, debug=debug)
}
else
{
  neglnl=sum((sapply(1:length(years), function(i,years){
     year=years[i]
     fS=gammad(par[(2*(i-1)+1):(2*i)],nmax)
     neglnl=-sum(log(prob[as.character(data$Start.year[seq(1,n,2)])==year,]%*%fS))
     if(!is.null(primary.only))
        neglnl=neglnl+ lnl.pods(par[(2*(i-1)+1):(2*i)], formula=formula, dpar=beta, 
            data=primary.only[as.character(primary.only$Start.year)==year,], gsS=gsS, debug=debug)
     return(neglnl)
     }, years=years)))
}
# if debug, output iteration results
if(debug)
{
   cat("\npar=",paste("c(",paste(par,collapse=","),")"))
   cat("\nneglnl=",neglnl)
}
return(neglnl)
}

fit.missed.pods <-
function(data, primary, pbyyear=FALSE, formula=~1, par=NULL, gsS,
                         maxit=1000,refit=TRUE,debug=FALSE,hessian=FALSE,method="BFGS")
{
# Fits true pod size distribution and the detection probability parameters
# for a given detection model specified by formula.
#
# Arguments:
#
#  data    - match data ordered in pairs with station=P then S where
#            P=Primary & S=Secondary
#  primary - observations made during primary period
#  pbyyear - if TRUE, fits year-specific gamma parameters for pod size
#  formula - the formula for the logistic detection model
#  par     - initial parameter values
#  gsS     - pod size calibration matrix
#  maxit   - maximum iterations
#  refit   - if TRUE will continue to call optim to refit until convergence is achieved
#             regardless of value of maxit
#  debug   - if TRUE will show iteration values
#  hessian - if TRUE will return hessian
#
#  Value:
#
#   list with elements : par  - parameter estimates
#                        AIC  - AIC value for model
#                        model- output from final optim run
################################################################################
nmax=ncol(gsS)
primary.only=primary[primary$only,]
if(nrow(primary.only)==0)primary.only=NULL
# add a dummy True pod size
data$True=1
# set number of pod size parameters depending on pbyyear and number of survey years
if(pbyyear)
{
   if(is.null(primary.only))
     years=levels(data$Start.year)
   else
     years=sort(unique(c(levels(primary.only$Start.year),levels(data$Start.year))))
   np=2*length(years)
}
else
{
   years=NULL
   np=2
}
# set number of parameters in detection probability model and create formula
# using podsize in place of True for io.glm
nc=ncol(model.matrix(formula,data))
if(length(grep("True",as.character(formula)))!=0)
   dformula=as.formula(paste("seen",paste(sub("True","podsize",as.character(formula)),collapse=""),sep=""))
else
   dformula=as.formula(paste("seen",paste(as.character(formula),collapse=""),sep=""))
# Depending on value of par argument create initial values for parameters
# Unless all are specified, a glm is used to specify the starting values for
# the detection parameters.  That is done using observed pod size (podsize) values
# rather than unknown true pod sizes.
if(is.null(par) | length(par)==np)
{
  if(length(par)==np)
    par=c(par,as.vector(coef(io.glm(data,dformula))))
  else
  {
    dpar=as.vector(coef(io.glm(data,dformula)))
    if(pbyyear)
    {
      psmod=podsize.computations(primary,gsS=gsS)
      spar=as.vector(sapply(psmod,function(x)x$par))
      sspar=NULL
      for(i in 1:length(years))
        sspar=c(sspar,fit.pods(spar[(2*(i-1)+1):(2*i)], formula=formula, dpar=dpar, data=primary[primary$Start.year==years[i],],  gsS=gsS, debug=FALSE, hessian=FALSE)$par)
    }
    else
    {
       sspar=c(.4,-.4)
    }
    par=c(sspar,dpar) 
  }
} else
{
  if(length(par)>nc+np)
    stop("initial par vector too long")
  else
    if(length(par) <nc+np) stop("initial vector too short")
}
if(debug) cat("\nInitial values = ",par,"\n")
# Fit model with optim which use lnl.missed.pods for the negative log-likelihood
mod=optim(par,lnl.missed.pods,years=years,formula=formula,data=data,primary.only=primary.only,
              gsS=gsS,control=list(maxit=maxit,ndeps=rep(1e-5,length(par))),debug=debug,hessian=hessian,method=method)
# If the model didn't converge and refit is TRUE, continue to fit the model
# using final values from last fit until it converges.
if(refit)
   while(mod$convergence!=0)
      mod=optim(mod$par,lnl.missed.pods,years=years,formula=formula,data=data,primary.only=primary.only,
                            gsS=gsS,control=list(maxit=maxit,ndeps=rep(1e-5,length(par))),debug=debug,hessian=hessian,method=method)
# Add names to each parameter in the model and return the result list
if(pbyyear)
{
  cnames=paste(rep(years,each=2),rep(c("Gamma shape","Gamma rate"),length(years)),sep=":")
  cnames=c(cnames,colnames(model.matrix(dformula,data)))
} else
  cnames=c("Gamma shape","Gamma rate",colnames(model.matrix(dformula,data)))
parvec=mod$par
names(parvec)=cnames
AIC=mod$value*2+2*length(cnames)
return(list(par=parvec,AIC=AIC,model=mod))
}

