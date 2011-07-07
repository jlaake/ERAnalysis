compute.series=function(models,naive.abundance,sightings=NULL,effort=NULL,
          gsS=NULL,best=TRUE,Match=NULL,cutoff=4,final.time=c(rep(90,20),100,100,90),
          lower.time=rep(0,23),lcut=0.2,mcut=1.0,twt=0.18,dwt=3.95,pwt=0.05,crittype="ellipse",
          DistBreaks=c(0,1,2,3,4,20),Use=TRUE,hessian=FALSE,debug=FALSE,TruePS=TRUE,
          recent.years=c(1987,1992,1993,1995,1997,2000,2001,2006),fn=1.0817)
{
# Define function to select best detection model and return list of models
# within a cutoff delta AICc value.
select.detection.models=function(x,models,cutoff=4)
{
   mod=vector("list",length(models))      
   for(i in 1:length(models))
      mod[[i]]=io.glm(x,as.formula(paste("seen~",models[i],sep="")))
   AICValues=sapply(mod,function(x) AIC(x))
   DeltaAIC=AICValues-min(AICValues)
   mod.numbers=(1:length(models))[order(DeltaAIC)]
   return(mod[mod.numbers[sort(DeltaAIC)<cutoff]])
}
# Get sightings and effort
# Assign day value and match up to Use field if Use=TRUE and select only those
# sightings with Use=TRUE.
PrimarySightings=NULL
PrimaryEffort=NULL
if(is.null(sightings))
{
   data(PrimarySightings,package="ERAnalysis",envir=environment())
   sightings=PrimarySightings
   PrimarySightings=NULL
}
if(is.null(effort))
{
   data(PrimaryEffort,package="ERAnalysis",envir=environment())
   effort=PrimaryEffort
   PrimaryEffort=NULL
}
sightings$day=as.Date(sightings$Date)-as.Date(paste(sightings$Start.year,"-12-01",sep=""))
sightings$seq=1:nrow(sightings)
if(Use)
{
   sightings=merge(sightings,subset(effort,select=c("key","Use")),by="key")
   sightings=sightings[sightings$Use,]
   sightings=sightings[order(sightings$seq),]
}
# if lcut >0, link sightings
if(lcut>0)
sightings=do.call("rbind",lapply(split(sightings, paste(sightings$Start.year,
        sightings$day, sep = "_")),
         function(x) return(linkdf(x,lcut=lcut,twt=twt,dwt=dwt,crittype=crittype))))
# Next unless the dataframe was provided create the Match dataframe with the 
#   specified parameters
if(is.null(Match))
   Match=create.match(lcut=lcut,mcut=mcut,twt=twt,dwt=dwt,pwt=pwt,crittype=crittype)
# Link Match to effort to obtain Use values if Use==TRUE
# in some cases the watch differs between primary and secondary matched detections
# when the t241 was at a boundary so Use value determined by Primary; also, in some
# cases there was no effort for primary because it was removed because
# there was no effort for the watch that met the beaufort/vis conditions.  In that
# case there is no matching record in primary effort so Use is set to FALSE.
if(Use)
{
   effort$key1=paste(effort$Date,effort$watch,sep="_")
   UseWatch=as.data.frame(table(effort$key1,effort$Use))
   UseWatch=UseWatch[UseWatch$Freq>0,1:2]
   names(UseWatch)=c("key","Use")
   Match=merge(Match,UseWatch,all.x=TRUE)
   Match$Use=factor(Match$Use,levels=c(TRUE,FALSE))
   Match$Use[is.na(Match$Use)]=FALSE
   Match=Match[order(Match$seq),]
   Match$Use=rep(Match$Use[Match$station=="P"],each=2)
   Match=Match[Match$Use==TRUE,]
}
# Run through the set of models for each of the recent years and select the
# best set of models within a delta aic "cutoff"
initial.models=vector("list",length(recent.years))
i=0
for (year in recent.years)
{
   i=i+1
   zz=Match[Match$Start.year==year,]
   zz$Start.year=factor(zz$Start.year)
   zz$Dist=cut(zz$distance,DistBreaks)
   zz$Observer=factor(zz$Observer)
   zz$Vis=cut(zz$vis,c(0,3,6))
   initial.models[[i]]=select.detection.models(zz,models,cutoff)
   print(summary(initial.models[[i]][[1]]))
   cat("\n",length(initial.models[[i]]))
}
# Run through the selected set of models and exclude any in which vis or
# beaufort effects are positive.
formulae=vector("list",length(recent.years))
for (i in 1:length(recent.years))
{
   for (j in 1:length(initial.models[[i]]))
   {
     est=coef(initial.models[[i]][[j]])
     vis.pos=grep("vis",names(est))
     Vis.pos=grep("Vis",names(est))
     beauf.pos=grep("beaufort",names(est))
     if(length(vis.pos)==0 &length(Vis.pos)==0 &length(beauf.pos)==0 )
       formulae[[i]]=c(formulae[[i]],formula(initial.models[[i]][[j]]))
     else
     {
        if(length(vis.pos)!=0 && est[vis.pos]<0)
          formulae[[i]]=c(formulae[[i]],formula(initial.models[[i]][[j]]))
        if(length(Vis.pos)!=0 && est[Vis.pos]<0)
          formulae[[i]]=c(formulae[[i]],formula(initial.models[[i]][[j]]))
        if(length(beauf.pos)!=0 && est[beauf.pos]<0)
          formulae[[i]]=c(formulae[[i]],formula(initial.models[[i]][[j]]))
     }
   }
}
# For each of the 8 recent years, compute the abundance estimate for each of the
# selected detection models unless best==TRUE
detection.models=vector("list",8)
abundance.models=vector("list",8)
#   Exclude effort and sightings for MAS in 1995 because she had no match data for 1995
sightings=sightings[!(sightings$Observer=="MAS"&sightings$Start.year==1995),]
effort=effort[!(effort$Observer=="MAS"&effort$Start.year==1995),]
i=0
if(is.null(gsS))data(gsS,package="ERAnalysis",envir=environment())
for (year in recent.years)
{
   i=i+1
#  Get match data and primary sightings for the year and set up variables
   zz=Match[Match$Start.year==year,]
   zz$Start.year=factor(zz$Start.year)
   primary=sightings[sightings$Start.year==year,]
   primary$Start.year=factor(primary$Start.year)
   primary$Dist=cut(primary$distance,DistBreaks)
   primary$Vis=cut(primary$vis,c(0,3,6))
   zz$Dist=cut(zz$distance,DistBreaks)
   zz$Vis=cut(zz$vis,c(0,3,6))
   zz$Observer=factor(zz$Observer)
   primary$Observer=factor(primary$Observer,levels=levels(zz$Observer))
#  Get effort data depending ion value of Use
   Start.year=NULL
   if(Use)
   {
      ern=subset(effort,subset=Use & as.character(Start.year)==year,
          select=c("Start.year","key","begin","end","effort","time","vis","beaufort"))
   } else
   {
      ern=subset(effort, as.character(Start.year)==year,
          select=c("Start.year","key","begin","end","effort","time","vis","beaufort"))
   }
   ern$Start.year=factor(ern$Start.year)
   if(i==1)
   {
      if(is.null(primary$corrected.podsize))   
         summary.df=data.frame(Year=year,nMatch=nrow(zz)/2,nPrimary=nrow(primary),
                            Whales=sum(primary$podsize),MeanPS=mean(primary$podsize),
                            HrsEffort=sum(ern$effort)*24)
      else
         summary.df=data.frame(Year=year,nMatch=nrow(zz)/2,nPrimary=nrow(primary),
                            Whales=sum(primary$podsize),MeanPS=mean(primary$podsize),CMeanPS=mean(primary$corrected.podsize),
                            HrsEffort=sum(ern$effort)*24)
   }
   else
   { 
      if(is.null(primary$corrected.podsize))   
         summary.df=rbind(summary.df,data.frame(Year=year,nMatch=nrow(zz)/2,nPrimary=nrow(primary),
                            Whales=sum(primary$podsize),MeanPS=mean(primary$podsize),
                            HrsEffort=sum(ern$effort)*24))
      else
         summary.df=rbind(summary.df,data.frame(Year=year,nMatch=nrow(zz)/2,nPrimary=nrow(primary),
                            Whales=sum(primary$podsize),MeanPS=mean(primary$podsize),CMeanPS=mean(primary$corrected.podsize),
                            HrsEffort=sum(ern$effort)*24))
   }
   cat("\n\nSurvey year               : ", year)
   cat("\n# Match records           : ", nrow(zz)/2)   
   cat("\n# Primary sightings       : ", nrow(primary))   
   cat("\n# Total whales            : ", sum(primary$podsize))   
   cat("\nMean observed pod size     : ", mean(primary$podsize))   
   if(!is.null(primary$corrected.podsize))   
   cat("\nMean corrected pod size    : ", mean(primary$corrected.podsize))   
   cat("\nHours of effort           : ", sum(ern$effort)*24)   
#  Loop over each detection model or just the best model
   nmodels=length(formulae[[i]])
   if(best) nmodels=1
   abundance.models[[i]]=vector("list",nmodels)
   detection.models[[i]]=vector("list",nmodels)
   for(j in 1:nmodels)
   {
     dformula=as.formula(paste("~", as.character(formulae[[i]][[j]])[3],sep=""))
#    If TruePS then replace podsize in formula with True which represents the unknown 
#      true podsize which is handled specifically by the code.
     if(TruePS & length(grep("podsize",as.character(dformula)))!=0)
        dtformula=as.formula(paste(sub("podsize","True",as.character(dformula)),collapse=""))
     else
        dtformula=dformula
#    Fit the detection/pod size model
     detection.models[[i]][[j]]=fit.missed.pods(formula=dtformula,pbyyear=TRUE,debug=debug,hessian=hessian,
          data=zz,primary=primary,gsS=gsS)
#    Compute the abundance estimate for this year; depends on whether TruePS=TRUE and
#    whether pod size was corrected to allow use of Reilly approach
     if(is.null(primary$corrected.podsize) | TruePS)
     {
       abundance.models[[i]][[j]]=estimate.abundance(spar=detection.models[[i]][[j]]$par[1:2],
         dpar=detection.models[[i]][[j]]$par[3:length(detection.models[[i]][[j]]$par)],gsS=gsS,effort=ern,
         sightings=primary, final.time=final.time[15+i],lower.time=lower.time[15+i],
         gformula=~s(time),dformula=dtformula)
     }
     else
     {
       abundance.models[[i]][[j]]=estimate.abundance(spar=NULL,
         dpar=detection.models[[i]][[j]]$par[3:length(detection.models[[i]][[j]]$par)],gsS=gsS,effort=ern,
         sightings=primary, final.time=final.time[15+i],lower.time=lower.time[15+i],
         gformula=~s(time),dformula=dtformula)
     }
     cat("\nEstimated whale abundance w/o fn: ",abundance.models[[i]][[j]]$Total)
   }
}
# Compute ratio of naive and corrected abundances for 1987-2006 (eqn 23)
Nbest=sapply(abundance.models,function(x) x[[1]]$Total)
Nnaive=sapply(naive.abundance[16:23],function(x)x$Total)
ratio=sum(Nbest)/sum(Nnaive)
# Compute series of estimates for 1967-2006 without nighttime correction factor (eqn 24)  
Nhat=c(sapply(naive.abundance[1:15],function(x)x$Total)*ratio, sapply(abundance.models,function(x) x[[1]]$Total))
# Apply nighttime correction factor (eqn 29)
Nhat=fn*Nhat
summary.df$Nhat=Nhat[16:23]
return(list(summary.df=summary.df,Match=Match,formulae=formulae,detection.models=detection.models,abundance.models=abundance.models,
                    ratio=ratio,Nhat=Nhat,TruePS=TruePS))
}

   