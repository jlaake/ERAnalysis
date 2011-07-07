compute.series.var=function(results,naive.abundance,ps.results,Use=TRUE,lcut=0.2,twt=0.18,dwt=3.95,crittype="ellipse",
                            final.time=c(rep(90,20),100,100,90), lower.time=rep(0,23),
                            DistBreaks=c(0,1,2,3,4,20),recent.years=c(1987,1992,1993,1995,1997,2000,2001,2006),
                            fn=1.0817,se.fn=0.0338,debug=FALSE,delta=c(0.001,0.01))
{
# Get PrimarySightings and effort
data(PrimarySightings,package="ERAnalysis",envir=environment())
data(PrimaryEffort,package="ERAnalysis",envir=environment())
data(gsS,package="ERAnalysis",envir=environment())
# Assign day value and match up to Use field if Use=TRUE
PrimarySightings$day=as.Date(PrimarySightings$Date)-as.Date(paste(PrimarySightings$Start.year,"-12-01",sep=""))
PrimarySightings$seq=1:nrow(PrimarySightings)
if(Use)
{
   PrimarySightings=merge(PrimarySightings,subset(PrimaryEffort,select=c("key","Use")),by="key")
   PrimarySightings=PrimarySightings[PrimarySightings$Use,]
   PrimarySightings=PrimarySightings[order(PrimarySightings$seq),]
}
# if lcut >0, link PrimarySightings
if(lcut>0)
PrimarySightings=do.call("rbind",lapply(split(PrimarySightings, paste(PrimarySightings$Start.year,
        PrimarySightings$day, sep = "_")),
         function(x) return(linkdf(x,lcut=lcut,twt=twt,dwt=dwt,crittype=crittype))))
#   Exclude effort and sightings for MAS in 1885 because she had no match data for 1995
PrimarySightings=PrimarySightings[!(PrimarySightings$Observer=="MAS"&PrimarySightings$Start.year==1995),]
PrimaryEffort=PrimaryEffort[!(PrimaryEffort$Observer=="MAS"&PrimaryEffort$Start.year==1995),]
#
# For each of the 8 recent years, compute the delta method variance due to variation in abundance estimate for each of the
# selected detection models unless best==TRUE
cat("Computing var1 component\n")
npar=0
for (i in 1:length(recent.years))
  npar=npar+ length(results$detection.models[[i]][[1]]$par)
vc.theta=matrix(0,nrow=npar,ncol=npar)
partial=matrix(0,nrow=23,ncol=npar)
i=0
jpos=0
Ns=sapply(results$abundance.models,function(x) x[[1]]$Total)
NaiveNs.late=sapply(naive.abundance[16:23],function(x)x$Total)
NaiveNs.early=sapply(naive.abundance[1:15],function(x)x$Total)
for (year in recent.years)
{
   i=i+1
#  Get effort data depending on value of Use
   Start.year=NULL
   if(Use)
   {
      ern=subset(PrimaryEffort,subset=Use & as.character(Start.year)==year,
          select=c("Start.year","key","begin","end","effort","time","vis","beaufort"))
   } else
   {
      ern=subset(PrimaryEffort, as.character(Start.year)==year,
          select=c("Start.year","key","begin","end","effort","time","vis","beaufort"))
   }
   ern$Start.year=factor(ern$Start.year)
   primary=PrimarySightings[PrimarySightings$Start.year==year,]
   primary$Start.year=factor(primary$Start.year)
   primary$Dist=cut(primary$distance,DistBreaks)
   primary$Vis=cut(primary$vis,c(0,3,6))
   primary$Observer=factor(primary$Observer,levels(factor(results$Match$Observer[results$Match$Start.year==year])))
   primary$seen=1
#  Loop over each parameter for group size and detection model and 
#  compute first derivatives of abundance estimates
   par.est=results$detection.models[[i]][[1]]$par
   vc.theta[(jpos+1):(jpos+length(par.est)),(jpos+1):(jpos+length(par.est))]=solve(results$detection.models[[i]][[1]]$model$hessian)
   dformula=as.formula(paste("~", as.character(results$formulae[[i]][[1]])[3],sep=""))
   if(results$TruePS & length(grep("podsize",as.character(dformula)))!=0)
      dtformula=as.formula(paste(sub("podsize","True",as.character(dformula)),collapse=""))
   else
      dtformula=dformula
   for (j in 1:length(par.est))
   {
     par.values=par.est
     par.values[j]=par.est[j]*(1-delta[1])
     nlower=estimate.abundance(spar=par.values[1:2],
       dpar=par.values[3:length(par.values)],gsS=gsS,effort=ern,
       sightings=primary, final.time=final.time[15+i],lower.time=lower.time[15+i],
       gformula=~s(time),dformula=dtformula,plotit=FALSE)$Total
     par.values[j]=par.est[j]*(1+delta[1])
     nupper=estimate.abundance(spar=par.values[1:2],
       dpar=par.values[3:length(par.values)],gsS=gsS,effort=ern,
       sightings=primary, final.time=final.time[15+i],lower.time=lower.time[15+i],
       gformula=~s(time),dformula=dtformula,plotit=FALSE)$Total
     partial[i+15,j+jpos]=(nupper-nlower)/(2*delta[1]*par.est[j])
     NewNs=Ns
     NewNs[i]=nlower
     ratio=sum(NewNs)/sum(NaiveNs.late)
     earlyNs.lower=ratio*NaiveNs.early
     NewNs[i]=nupper
     ratio=sum(NewNs)/sum(NaiveNs.late)
     earlyNs.upper=ratio*NaiveNs.early
     partial[1:15,j+jpos]=(earlyNs.upper-earlyNs.lower)/(2*delta[1]*par.est[j])
   }
   jpos=jpos+length(par.est)
}
# This is the first part of the v-c matrix for the N series
var1=partial%*%vc.theta%*%t(partial)
cat("Completed var1 component\n")
# Next compute the v-c matrix due to variation in pod size calibration data
par.est=ps.results$par
vc.theta=ps.results$vc
partial=matrix(0,nrow=23,ncol=length(par.est))
cat("Computing var2 component\n")
for (j in 1:length(par.est))
{
cat("Parameter ",j,"\n")
#    par-delta*par
  par.values=par.est
  par.values[j]=par.est[j]*(1-delta[2])
  nmax=20
  ps.results$par=par.values
  gsS=create.gsS(ps.results,nmax=20)  
  Nlower=rep(0,23)
  i=0
  for (year in recent.years)
  {
    i=i+1
#   Get effort data depending on value of Use
    if(Use)
    {
      ern=subset(PrimaryEffort,subset=Use & as.character(Start.year)==year,
          select=c("Start.year","key","begin","end","effort","time","vis","beaufort"))
    } else
    {
      ern=subset(PrimaryEffort, as.character(Start.year)==year,
          select=c("Start.year","key","begin","end","effort","time","vis","beaufort"))
    }
    ern$Start.year=factor(ern$Start.year)
    zz=results$Match[results$Match$Start.year==year,]
    zz$Start.year=factor(zz$Start.year)
    zz$Dist=cut(zz$distance,DistBreaks)
    zz$Vis=cut(zz$vis,c(0,3,6))
    zz$Observer=factor(zz$Observer)
    primary=PrimarySightings[PrimarySightings$Start.year==year,]
    primary$Start.year=factor(primary$Start.year)
    primary$Dist=cut(primary$distance,DistBreaks)
    primary$Vis=cut(primary$vis,c(0,3,6))
    primary$Observer=factor(primary$Observer,levels(factor(zz$Observer)))
    primary$seen=1
    dformula=as.formula(paste("~", as.character(results$formulae[[i]][[1]])[3],sep=""))
    if(results$TruePS & length(grep("podsize",as.character(dformula)))!=0)
      dtformula=as.formula(paste(sub("podsize","True",as.character(dformula)),collapse=""))
    else
      dtformula=dformula
    dpar=results$detection.models[[i]][[1]]$par
    ddpar=fit.missed.pods(formula=dtformula,pbyyear=TRUE,debug=FALSE,hessian=FALSE,par=dpar,
          data=zz,primary=primary,gsS=gsS)
    if(debug)
    {
       cat("\nconvergence =",ddpar$model$convergence)
       cat("\nvalue =",ddpar$model$value)
       cat("\npar =",ddpar$model$par)
       cat("\ncount =",ddpar$model$counts)
    }
    ddpar=ddpar$par
    spar=ddpar[1:2]
    ddpar=ddpar[3:length(ddpar)]
    Nlower[i+15]=estimate.abundance(spar=spar,dpar=ddpar,gsS=gsS,effort=ern,
       sightings=primary, final.time=final.time[15+i],lower.time=lower.time[15+i],
       gformula=~s(time),dformula=dtformula,plotit=FALSE)$Total
  }
    if(debug)cat("\nNlower=",Nlower)
    par.values=par.est
    par.values[j]=par.est[j]*(1+delta[2])
    ps.results$par=par.values
    gsS=create.gsS(ps.results,nmax=20)  
    Nupper=rep(0,23)
    i=0
    for (year in recent.years)
    {
      i=i+1
#     Get effort data depending on value of Use
      if(Use)
      {
        ern=subset(PrimaryEffort,subset=Use & as.character(Start.year)==year,
          select=c("Start.year","key","begin","end","effort","time","vis","beaufort"))
      } else
      {
        ern=subset(PrimaryEffort, as.character(Start.year)==year,
          select=c("Start.year","key","begin","end","effort","time","vis","beaufort"))
      }
      ern$Start.year=factor(ern$Start.year)
      zz=results$Match[results$Match$Start.year==year,]
      zz$Start.year=factor(zz$Start.year)
      zz$Dist=cut(zz$distance,DistBreaks)
      zz$Vis=cut(zz$vis,c(0,3,6))
      zz$Observer=factor(zz$Observer)
      primary=PrimarySightings[PrimarySightings$Start.year==year,]
      primary$Start.year=factor(primary$Start.year)
      primary$Dist=cut(primary$distance,DistBreaks)
      primary$Vis=cut(primary$vis,c(0,3,6))
      primary$Observer=factor(primary$Observer,levels(factor(zz$Observer)))
      primary$seen=1
      dformula=as.formula(paste("~", as.character(results$formulae[[i]][[1]])[3],sep=""))
      if(results$TruePS & length(grep("podsize",as.character(dformula)))!=0)
        dtformula=as.formula(paste(sub("podsize","True",as.character(dformula)),collapse=""))
      else
        dtformula=dformula
      dpar=results$detection.models[[i]][[1]]$par
      ddpar=fit.missed.pods(formula=dtformula,pbyyear=TRUE,debug=FALSE,hessian=FALSE,par=dpar,
          data=zz,primary=primary,gsS=gsS)
      if(debug)
      {
         cat("\nconvergence =",ddpar$model$convergence)
         cat("\nvalue =",ddpar$model$value)
         cat("\npar =",ddpar$model$par)
         cat("\ncount =",ddpar$model$counts)
      }
      ddpar=ddpar$par
      spar=ddpar[1:2]
      ddpar=ddpar[3:length(ddpar)]
      Nupper[i+15]=estimate.abundance(spar=spar,dpar=ddpar,gsS=gsS,effort=ern,
         sightings=primary, final.time=final.time[15+i],lower.time=lower.time[15+i],
         gformula=~s(time),dformula=dtformula,plotit=FALSE)$Total
    }
    if(debug) cat("\nNupper=",Nupper)
    partial[16:23,j]=(Nupper[16:23]-Nlower[16:23])/(2*delta[2]*par.est[j])
    ratio=sum(Nlower[16:23])/sum(NaiveNs.late)
    earlyNs.lower=ratio*NaiveNs.early
    ratio=sum(Nupper[16:23])/sum(NaiveNs.late)
    earlyNs.upper=ratio*NaiveNs.early
    partial[1:15,j]=(earlyNs.upper-earlyNs.lower)/(2*delta[2]*par.est[j])
  }
if(debug)write.table(partial,"partial2.txt")
var2=partial%*%vc.theta%*%t(partial)
cat("Completed var2 component\n")
# V-c matrix component from fitting gam and residual variance around gam
cat("Computing var3 component\n")
var3=matrix(0,nrow=23,ncol=23)
var.Nhat=sapply(results$abundance.models,function(x) x[[1]]$var.Total)
diag(var3[16:23,16:23])=var.Nhat
ratio=sum(Ns)/sum(NaiveNs.late)
sigma.ratio=(var(Ns)+ratio^2*var(NaiveNs.late)-2*ratio*cov(Ns,NaiveNs.late))/(length(recent.years)*mean(NaiveNs.late)^2)
var.NaiveNs=sapply(naive.abundance, function(x) x$var.Total)
var.NaiveNs.early=var.NaiveNs[1:15]
var.NaiveNs.late=var.NaiveNs[16:23]
var3[1:15,1:15]=sigma.ratio*outer(NaiveNs.early,NaiveNs.early,"*")
diag(var3[1:15,1:15])=ratio^2*NaiveNs.early^2*
              (sigma.ratio*(1+length(recent.years))/ratio^2+var.NaiveNs.early/NaiveNs.early^2)
covar.ratio=outer(NaiveNs.early,(var.Nhat/NaiveNs.late-Ns*ratio*var.NaiveNs.late/NaiveNs.late^2),"*")
var3[1:15,16:23]=covar.ratio
var3[16:23,1:15]=t(covar.ratio)
vc=var1+var2+var3
se=sqrt(diag(vc))
# Add adjustment for nighttime correction factor; first line does covariances f^2*cov(Nhat_i,Nhat_j)
vc=vc+fn^2*vc
# Following adjusts variances
diag(vc)=(results$Nhat)^2*((se.fn/fn)^2+(se/results$Nhat)^2)
se=sqrt(diag(vc))
return(list(vc=vc,var1=var1,var2=var2,var3=var3,se=se,cormat=vc/outer(se,se,"*")))
}
