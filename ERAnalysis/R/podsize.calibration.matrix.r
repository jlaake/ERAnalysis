create.podsize.calibration.matrix=function(package.dir="C:/Users/Jeff Laake/Desktop/MyDocuments/R Development/ERAnalysis/ERAnalysis/data",save=TRUE)
{
################################################################################
# Estimate pod size calibration matrix (gsS) from pod size calibration data and fitted
# gamma.pod model (assumed to be in workspace) and create additive 
# correction factor vectors and store all as data for package
################################################################################
PodsizeCalibrationData=NULL
PodsizeCalibrationTable=NULL
data(PodsizeCalibrationData,package="ERAnalysis",envir=environment())
data(PodsizeCalibrationTable,package="ERAnalysis",envir=environment())
# Create Reilly additive correction factors
# add.cf.all pools all of the data
# add.cf.reilly only uses 1970s aerial data that was applied in Buckland et al 1992
# add.cf.laake only uses 1990s aerial data that was applied in all surveys from 1992 forward
nmax=20
psdf=PodsizeCalibrationTable
add.cf.all=with(PodsizeCalibrationData,tapply(True-Estimate,cut(Estimate,breaks=c(0,1,2,3,20)),mean))
add.cf.reilly=with(PodsizeCalibrationData[PodsizeCalibrationData$Year=="1978/1979",],tapply(True-Estimate,cut(Estimate,breaks=c(0,1,2,3,20)),mean))
add.cf.laake=with(PodsizeCalibrationData[PodsizeCalibrationData$Type=="Aerial"&PodsizeCalibrationData$Year!="1978/1979",],
       tapply(True-Estimate,cut(Estimate,breaks=c(0,1,2,3,20)),mean))
if(save)save(add.cf.all,file=file.path(package.dir,"add.cf.all.rda"))
if(save)save(add.cf.reilly,file=file.path(package.dir,"add.cf.reilly.rda"))
if(save)save(add.cf.laake,file=file.path(package.dir,"add.cf.laake.rda"))
# Create gsS matrix from pod random effects model with gamma distribution
gss=create.gsS(gamma.pod,nmax=20)
row.names(gsS)=NULL
if(save)save(gsS,file=file.path(package.dir,"gsS.rda"))
return(NULL)
}

create.gsS=function(ps.results,nmax=20,True=NULL)
{
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
    return(dnorm(x,sd=sigma)*gam.d(shape,rate,x=rep(obs,length(x)),nmax=nmax))
  }
  gpar=ps.results$par
  sigma=exp(gpar[1])
  if(is.null(True))
    nprime=nmax
  else
    nprime=length(True)
  gsS=matrix(0,nrow=nprime,ncol=nmax)
  for (i in 1:nprime)
  {
     if(is.null(True))
        x=data.frame(size=cut(i,c(1,2,3,4,nmax+1),right=FALSE),plus=as.numeric(i>3),True=i)
     else
        x=data.frame(size=cut(True[i],c(1,2,3,4,nmax+1),right=FALSE),plus=as.numeric(True[i]>3),True=True[i])     
     shape=exp(model.matrix(~size,x)%*%gpar[7:10])
     xbeta=model.matrix(~size+True:plus,x)%*%gpar[2:6]
     for (j in 1:nmax)
       gsS[i,j]=integrate(gam.reint,-5*sigma,5*sigma,xbeta=xbeta,shape=shape,sigma=sigma,obs=j)$val
  }
  return(gsS)
}
################################################################################
reilly.cf=function(x,add.cf) x+add.cf[as.numeric(cut(x,breaks=c(0,1,2,3,20)))]
################################################################################

################################################################################
expected.podsize=function(spar,gsS,nmax=20)
{
# Compute true pod size distribution and ES
  fS=gammad(spar,nmax)
  ES=sum(fS*(1:nmax))
# Compute conditional distribution and conditional expectation for
# true pod size given an observed pod size
  fSs=fS*gsS/matrix(colSums(fS*gsS),nrow=nmax,ncol=nmax,byrow=TRUE)
  ESs=apply(fSs,2,function(x) sum((1:nmax)*x))
  return(list(fS=fS,ES=ES,fSs=fSs,ESs=ESs))
}
################################################################################

podsize.computations=function(x,gsS,nmax=20,plot=FALSE)
{
################################################################################
# Observed pod size distribution and estimation of true pod size distribution
# assuming equal visibility of all sizes; returns parameters for "true" pod
# size distribution which can be used as starting values for the more
# involved computations of the parameters of true podsize.
################################################################################
# Create observed pod size table and summary of number of pods used
ps.table=table(x$Start.year,factor(x$podsize,levels=1:nmax))
# Plot combined pod size distribution based on proportions
all.years=unique(x$Start.year)
if(plot)
{
  dev.new()
  barplot(ps.table/rowSums(ps.table),beside=TRUE,col=heat.colors(nrow(ps.table)),
         legend.text=all.years)
# Plot observed pod size distributions
  pdf("ObservedPodsizeDistributions.pdf")
  par(mfrow=c(3,2))
  for(i in 1:length(all.years))
     barplot(ps.table[i,],xlab="Pod size",xlim=c(1,nmax),main=all.years[i])
  dev.off()
}
# For each year fit the distribution and plot the fitted/observed distribution
# This assumes that all sizes are equally visible which is not the case
# The fitted models are stored in the list psmod
if(plot)
{
  pdf("fittedps.pdf")
  par(mfrow=c(4,3),mar=c(2,2,1,2))
}
psmod=vector("list",length=length(all.years))
fS=matrix(0,nrow=length(all.years),ncol=nmax)
i=0
for(year in all.years)
{
  i=i+1
  psmod[[i]]=optim(c(0,0),fit.fS,gsS=gsS,pstable=ps.table[i,])
  fS[i,]=gammad(psmod[[i]]$par,nmax)
  if(plot)barplot(rbind(colSums(gsS*fS[i,])*sum(ps.table[i,]),ps.table[i,]),beside=TRUE,legend.text=c("Exp","Obs"),main=paste("Year=",year,"/",year+1,sep=""))
}
# Next plot the estimated "true" distribution of pod size and the distribution of
# measured (observed) pod size.
if(plot)
{
  dev.off()
  pdf("fS.pdf")
  par(mfrow=c(4,3),mar=c(2,2,1,2))
  i=0
  for(year in all.years)
  {
    i=i+1
    barplot(rbind(fS[i,]*sum(ps.table[i,]),ps.table[i,]),beside=TRUE,legend.text=c("True","Estimate"),main=paste("Year=",year,"/",year+1,sep=""))
  }
  dev.off()
}
return(psmod)
}

