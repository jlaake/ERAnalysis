estimate.abundance <-
function(spar,dpar,gsS,effort,sightings,dformula=~True,gformula=~s(time),nmax=20,
          pod=FALSE,plotit=TRUE,anchor=TRUE,show.anchor=FALSE,sp=NULL,final.time=90,
          lower.time=0,do.mult=TRUE,pool=TRUE,...)
{
# Function to compute observation-specific detection probability
# from formula parameters and sightings data
probs=function(j)
{
  sightings$True=j
  dmat=model.matrix(dformula,sightings)
  plogis(dmat%*%dpar)
}
# Order sightings by Start.year and create list of years
sightings=sightings[order(sightings$Start.year),]
years=levels(factor(sightings$Start.year))
Nhat.whales=NULL
Nhat=NULL
for(year in years)
{
   cat("\n\nSurvey year               : ", year)
   cat("\n# Primary sightings       : ", nrow(sightings[sightings$Start.year==year,]))   
   cat("\n# Total whales            : ", sum(sightings$podsize[sightings$Start.year==year]))   
   cat("\nMean observed podsize     : ", mean(sightings$podsize[sightings$Start.year==year]))   
   cat("\nHours of effort           : ", sum(effort$effort)*24)   
}
# Depending on values of dpar and spar construct estimates of
# number of whales and number of pods passing for each observation
if(!is.null(dpar)&!is.null(dformula)&!is.null(spar))
{
   # Compute matrix of detection probabilities for each observation (rows) at each
   # of 1 to nmax possible true pod sizes
   ps=sapply(1:nmax, probs)
   # Corrections for pod size and detection probability
   # Compute probability function for True pod size
   i=0                  
   for(year in years)
   {
      i=i+1
      fS=gammad(spar[((i-1)*2+1):(i*2)],nmax)
      # Compute conditional detection probability function for true size given observed size
      fSs=t(t(fS*gsS)/colSums(fS*gsS))
      # Compute estimate of expected number of whales represented by observed whales
      Nhat.whales=c(Nhat.whales,rowSums(t(fSs[,sightings$podsize[sightings$Start.year==year]]*(1:nmax))/ps[sightings$Start.year==year]))
      # Same as above for pods rather than whales
      Nhat=c(Nhat,rowSums(t(fSs[,sightings$podsize[sightings$Start.year==year]])/ps[sightings$Start.year==year]))
   }
} else
{
  # Corrections for detection probability but
  # no further corrections for pod size but pod size
  # may have been corrected prior as in reilly.cf
  if(!is.null(dpar)&!is.null(dformula))
  {
    dmat=model.matrix(dformula,sightings)
    ps=plogis(dmat%*%dpar)
    if(is.null(sightings$corrected.podsize))
       Nhat.whales=sightings$podsize/ps
      else
       Nhat.whales=sightings$corrected.podsize/ps
    Nhat=1/ps
  }else
  {
    # Corrections for pod size but no corrections
    # detection probability
    if(!is.null(spar))
    {
       i=0
      for(year in years)
      {
         i=i+1                                                       
         # Compute probability function for True pod size
         fS=gammad(spar[((i-1)*2+1):(i*2)],nmax)
         # Compute conditional detection probability function for true size given observed size
         fSs=t(t(fS*gsS)/colSums(fS*gsS))
         # Compute estimate of expected number of whales represented by observed whales
         Nhat.whales=c(Nhat.whales,rowSums(t(fSs[,sightings$podsize[sightings$Start.year==year]]*(1:nmax))))
      }
      Nhat=rep(1,nrow(sightings))
    } else
    # No further corrections for either but pod size
    # may have been corrected prior as in reilly.cf
    {
      if(is.null(sightings$corrected.podsize))
         Nhat.whales=sightings$podsize
      else
         Nhat.whales=sightings$corrected.podsize
      Nhat=rep(1,nrow(sightings))
    }
  }
}
#
# Merge estimates with effort data to construct migration timing model
# to extrapolate from sampled periods to entire migration time period
#
sightings$Nhat.whales=Nhat.whales
sightings$Nhat=Nhat
if(!pod)
  est.df=data.frame(key=unique(sightings$key),
          nhat=sapply(split(sightings$Nhat.whales,factor(sightings$key)),sum))
else
  est.df=data.frame(key=unique(sightings$key),
          nhat=sapply(split(sightings$Nhat,factor(sightings$key)),sum))
ern=subset(effort,select=c("Start.year","key","begin","end","effort","time","vis","beaufort"))
ern=merge(ern,est.df,by.x="key",by.y="key",all.x=TRUE)
ern$nhat[is.na(ern$nhat)]=0
results=fit.migration.gam(ern,years=as.numeric(years),formula=gformula,pod=pod,plotit=plotit,
                            anchor=anchor,show.anchor=show.anchor,sp=sp,
                            final.time=final.time,lower.time=lower.time,do.mult=do.mult,pool=pool,...)
results$options=list(anchor=anchor,pod=pod,final.time=final.time,lower.time=lower.time,
                     gformula=gformula,dformula=dformula,spar=spar,dpar=dpar,nmax=nmax)
cat("\nEstimated whale abundance : ",results$Total)
return(results)
}

