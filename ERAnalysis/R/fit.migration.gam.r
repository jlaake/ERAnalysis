fit.migration.gam=function(er.migdata,years,formula=~s(time),pod=FALSE,plotit=TRUE,anchor=TRUE,
                                show.anchor=FALSE,sp=NULL,final.time=rep(90,length(years)),
                                lower.time=rep(0,length(years)),do.mult=TRUE,pool=TRUE,...)
{
   final=final.time
   lower=lower.time
   if(length(final)!=length(years))stop("Number of elements in final.time does not match number of years")
   if(length(lower)!=length(years))stop("Number of elements in lower.time does not match number of years")
   require(mgcv)
   ermod=vector("list",length(years))
   formula=as.formula(paste("nhat",paste(as.character(formula),collapse=""),sep=""))
   if(pod)
     ylabel="Pods per day"
   else
     ylabel="Whales per day"
   if(!pool)
   {
   i=0
   for (year in years)
   {
     i=i+1
     er=er.migdata[er.migdata$Start.year==year,]
     er$offset=log(er$effort)
     if(anchor)
     {
        if(lower[i] <= (er$begin[1]-1))
        {
           er=rbind(er[1,],er)
           er$time[1]=c(lower[i]+0.5)
           er$begin[1]=c(lower[i])
           er$end[1]=c(lower[i]+1)
           er$offset[1]=c(0)
           er$effort[1]=1
           er$vis[1]=1
           er$beaufort[1]=0           
           er$nhat[1]=0
        }
        else
           lower[i]=floor(er$begin[1])
        if(er$end[nrow(er)]<(final[i]-1))
        {
           er=rbind(er,er[1,])
           er$time[nrow(er)]=final[i]-0.5
           er$begin[nrow(er)]=final[i]-1
           er$end[nrow(er)]=final[i]
           er$offset[nrow(er)]=0
           er$effort[nrow(er)]=1
           er$nhat[nrow(er)]=0
           er$vis[nrow(er)]=1
           er$beaufort[nrow(er)]=0           

        }
        else
           final[i]=ceiling(er$end[nrow(er)])
     }
     ermod[[i]]=gam(formula,data=er,offset=offset,family=quasipoisson)
     ppd=tapply(er$nhat,floor(er$time),sum)/tapply(er$effort,floor(er$time),sum)
     if(plotit)
     {
        Eppd=predict(ermod[[i]],type="response")
        ymax=max(c(Eppd,ppd))*1.05
        plot(er$time,Eppd, ylim=c(0,ymax),type="l",main=paste(year,"/",year+1,sep=""),xlab="Days since 1 Dec", ylab=ylabel,xlim=c(0,100))
        points(as.numeric(names(ppd)),ppd)
        if(anchor & show.anchor)
        {
           if(lower.time[i] <= (er$begin[1]-1)) er=er[-1,]
           if(er$end[nrow(er)-1]<(final.time[i]-1))er=er[-nrow(er),]
           mod=gam(formula,data=er,offset=offset,family=quasipoisson)
           newdata=data.frame(time=(floor(lower[i]):(final[i]-1))+.5)
           newdata$Start.year=factor(rep(year,nrow(newdata)),levels=years)
           lines((floor(lower[i]):(final[i]-1))+.5,
             predict(mod,newdata=newdata,type="response"),
             xlim=c(lower[i],final[i]),lty=2)
        }
     }
   }
   Total=vector("numeric",length(years))
   for(i in 1:length(years))
   {
        newdata=data.frame(time=(lower[i]:(final[i]-1))+.5)
        newdata$Start.year=factor(rep(years[i],nrow(newdata)),levels=years)
        if(length(grep("vis",as.character(formula)))!=0)
          newdata$vis=1
        if(length(grep("beaufort",as.character(formula)))!=0)
          newdata$beaufort=0
        Total[i]=sum(predict(ermod[[i]],newdata=newdata,type="response"))
   }
   if(do.mult)
   {
      TotalM=vector("numeric",length(years))
      mult=vector("numeric",length(years))
      i=0
      for(year in years)
      {
        i=i+1
        mult[i]=compute.sampling.multiplier(ermod[[i]],er.migdata[er.migdata$Start.year==year,],
                upper=final[i],lower=lower[i])
        TotalM[i]=sum(ermod[[i]]$y)*mult[i]
      }
      return(list(models=ermod,Total=Total,TotalM=TotalM,mult=mult))
   }
   else
      return(list(models=ermod,Total=Total))
   }
   else
   {
   i=0
   ern=NULL
   for (year in years)
   {
     i=i+1
     er=er.migdata[er.migdata$Start.year==year,]
     er$offset=log(er$effort)
     if(anchor)
     {
        if(lower[i] <= (er$begin[1]-1))
        {
           er=rbind(er[1,],er)
           er$time[1]=c(lower[i]+0.5)
           er$begin[1]=c(lower[i])
           er$end[1]=c(lower[i]+1)
           er$offset[1]=c(0)
           er$effort[1]=1
           er$nhat[1]=0
           er$vis[1]=1
           er$beaufort[1]=0                      
        }
        else
           lower[i]=floor(er$begin[1])
        if(er$end[nrow(er)]<(final[i]-1))
        {
           er=rbind(er,er[1,])
           er$time[nrow(er)]=final[i]-0.5
           er$begin[nrow(er)]=final[i]-1
           er$end[nrow(er)]=final[i]
           er$offset[nrow(er)]=0
           er$effort[nrow(er)]=1
           er$nhat[nrow(er)]=0
           er$vis[nrow(er)]=1
           er$beaufort[nrow(er)]=0           
        }
        else
           final[i]=ceiling(er$end[nrow(er)])
     }
     ern=rbind(ern,er)
     }
     ern$Start.year=factor(ern$Start.year)
     if(!is.null(sp))
     {
        ermod=gam(formula,data=ern,offset=offset,family=quasipoisson,sp=sp,...)
        ermod=gam(formula,data=ern,offset=offset,family=quasipoisson,start=ermod$coef,...)
     } else
        ermod=gam(formula,data=ern,offset=offset,family=quasipoisson,...)     
     if(anchor & show.anchor)
     {
       er=er.migdata
       er$offset=log(er$effort)
       er$Start.year=factor(er$Start.year)
       ermod.na=gam(formula,data=er,offset=offset,family=quasipoisson,...)
       Eppd.all.na=predict(ermod.na,newdata=ern,type="response")
     }
     i=0
     Eppd.all=predict(ermod,type="response")
     for (year in years)
     {
        i=i+1
        er=ern[ern$Start.year==year,]
        ppd=tapply(er$nhat,floor(er$time),sum)/tapply(er$effort,floor(er$time),sum)
        if(plotit)
        {
           Eppd=Eppd.all[ern$Start.year==year]
           ymax=max(c(Eppd,ppd))*1.05
           plot(er$time,Eppd, ylim=c(0,ymax),type="l",main=paste(year,"/",year+1,sep=""),xlab="Days since 1 Dec", ylab=ylabel,xlim=c(0,100))
           points(as.numeric(names(ppd)),ppd)
           if(anchor & show.anchor)
           {
           lines((floor(lower[i]):(final[i]-1))+.5,
             Eppd.all.na[ern$Start.year==year],
             xlim=c(lower[i],final[i]),lty=2)
        }
       }
     }
     newdata=NULL
     for(i in 1:length(years))
     {
        ndata=data.frame(time=(lower[i]:(final[i]-1))+.5)
        ndata$Start.year=factor(rep(years[i],nrow(ndata)),levels=years)
        if(length(grep("vis",as.character(formula)))!=0)
          ndata$vis=1
        if(length(grep("beaufort",as.character(formula)))!=0)
          ndata$beaufort=0
        newdata=rbind(newdata,ndata)
     }
     Eppd.all=predict(ermod,newdata=newdata,type="response")
     Total=tapply(Eppd.all,newdata$Start.year,sum)  
     var.Total=NULL 
     sumN=NULL 
     for(i in 1:length(years))
     {
       newd=newdata[newdata$Start.year==years[i],]
       Xp <- predict(ermod,newd,type="lpmatrix") 
       Xs=Xp*as.vector(exp(Xp%*%coef(ermod)))
       var.Total=c(var.Total,sum(Xs%*%ermod$Vp%*%t(Xs))+sum(Eppd.all*ermod$scale))
      }   
     return(list(models=ermod,Total=Total,var.Total=var.Total,pred=Eppd.all))
   }
}

        