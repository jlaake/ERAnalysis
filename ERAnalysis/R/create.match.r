create.match=function(lcut=0.2,mcut=1,twt=.18,dwt=3.95,pwt=0.05,crittype="ellipse",fdir="c:/gw")
{
  Ld2Ll  <- function(Ld) 
  {
     Ll = split(Ld, paste(Ld$Start.year,Ld$day, sep = "_"))
     return(Ll)
  }
# Extract match data and call algorithm to make links/matches
  Match.list<- extract.match.data()
  if(lcut>0)
  {
     Linked.data=do.call("rbind",lapply(Match.list,
         function(x) return(linklist(x,lcut=lcut,twt=twt,dwt=dwt,crittype=crittype))))
     Linked.list <- Ld2Ll(Linked.data)
     Matched.data=do.call("rbind",lapply(Linked.list,function(x) 
         return(matchdf(x,mcut=mcut,twt=twt,dwt=dwt,pwt=pwt,crittype=crittype))))
  }
  else
  {
     Matched.data=do.call("rbind",lapply(Match.list,function(x) 
         return(matchdf(x,mcut=mcut,twt=twt,dwt=dwt,pwt=pwt,crittype=crittype))))

  }
# Rename and select needed fields and create key field
  Matched.data$station="P"
  Matched.data$station[Matched.data$PorS==2]="S"
  Matched.data$station=factor(Matched.data$station) 
  Match=NULL
  AllPrimaryEffort=NULL
  AllSecondaryEffort=NULL
  data(Match,package="ERAnalysis",envir=environment())
  match.names=names(Match)
  xx=subset(Matched.data,select=c(names(Matched.data)[names(Matched.data)%in%match.names],"observer"))
  xx$Date=as.character(strptime(paste(xx$Start.year,"-12-01",sep=""),format="%Y-%m-%d")+3600*24*(xx$day-1))
  xx$key=paste(xx$Date,xx$watch,sep="_")
# get observer, beaufort, vis
# first do Primary
  data(AllPrimaryEffort,package="ERAnalysis",envir=environment())
  data(AllSecondaryEffort,package="ERAnalysis",envir=environment())
  PrimaryEffort=AllPrimaryEffort
  xx$seq=1:nrow(xx)
  zz=merge(xx[xx$station=="P",],PrimaryEffort,by="Date")
  zz$begin.time=24*(zz$begin-floor(zz$begin))
  zz$end.time=24*(zz$end-floor(zz$end))
  zz=zz[zz$t241<=(zz$end.time+1/3600)&zz$t241>=zz$begin.time,]
  zz=zz[order(zz$seq),]
# next do Secondary
  SecondaryEffort=AllSecondaryEffort
  zz.sec=merge(xx[xx$station=="S",],SecondaryEffort,by="Date")
  zz.sec$begin.time=24*(zz.sec$begin-floor(zz.sec$begin))
  zz.sec$end.time=24*(zz.sec$end-floor(zz.sec$end))
  zz.sec=zz.sec[zz.sec$t241<=(zz.sec$end.time+1/3600)&zz.sec$t241>=zz.sec$begin.time,]
  zz.sec=zz.sec[order(zz.sec$seq),]
# create list of records in which t241 of both sightings are in effort period
  on.effort=(xx$seq[xx$station=="P"] %in% zz$seq) & (xx$seq[xx$station=="S"] %in% zz.sec$seq)
  xx=xx[rep(on.effort,each=2),]
  zz=zz[zz$seq %in% xx$seq,]
  zz.sec=zz.sec[zz.sec$seq %in% xx$seq,]
# extract observer record from effort record containing t241;  this changes observer
# on record of about 2-3% of the sightings
  xx$observer[xx$station=="P"]=as.character(zz$Observer)
  xx$observer[xx$station=="S"]=as.character(zz.sec$Observer)
  xx$observer[xx$observer=="35"]="DJR"
# add beaufort and vis
  xx$beaufort[xx$station=="P"]=zz$beaufort.y
  xx$beaufort[xx$station=="S"]=zz.sec$beaufort.y
  xx$vis[xx$station=="P"]=zz$vis.y
  xx$vis[xx$station=="S"]=zz.sec$vis.y
# create observer experience table and add hours of experience to match data
  ObserverExp=CreateObserverExperience()
  xx=merge(xx,ObserverExp,by="Date",all.x=TRUE)
  xx$observer=as.character(xx$observer)
  xx$hours=xx[cbind(1:nrow(xx),match(xx$observer,names(xx)[16:74])+15)]
  xx=xx[,c(1:15,75)]
  xx$hours=as.numeric(xx$hours)
# next construct pods per hour (as determined from primary observer) during the watch
# containing the match data
  PrimaryEffort$key=paste(PrimaryEffort$Date,PrimaryEffort$watch,sep="_")
  pphr=with(PrimaryEffort, tapply(npods,key,sum))
  pphr=pphr/with(PrimaryEffort, tapply(effort,key,sum))
  xx=merge(xx,data.frame(key=names(pphr),pphr=pphr),by="key")
# Finally add on the sex of the observer
  data(Observer,package="ERAnalysis",envir=environment())
  Observer$key=as.character(Observer$Initials)
  Observer$key[is.na(Observer$Initials)]=Observer$Observer[is.na(Observer$Initials)]
  xx=merge(xx, subset(Observer, select = c("key", "Sex")), by.x ="observer",by.y = "key", all.x = TRUE)
  xx$Observer=factor(xx$observer)
  xx$observer=NULL
  xx=xx[order(xx$seq),]
  row.names(xx)=NULL
  return(xx)
}
CreateObserverExperience=function()
{
# Create dataframe of hours of experience for each observer
# First extract values from early and recent data
EarlyEffort=NULL
ERSurveyData=NULL
data(EarlyEffort,package="ERAnalysis",envir=environment())
EarlyExp=data.frame(Observer=EarlyEffort$Observer,Date=substr(EarlyEffort$key,1,10),Hours=as.vector(EarlyEffort$End.date.time-EarlyEffort$Begin.date.time))
data(ERSurveyData,package="ERAnalysis",envir=environment())
EFLAG=NULL
EXPERIMENT=NULL
Start.watch=subset(ERSurveyData,subset=EFLAG==1&EXPERIMENT%in%c(1,2),select=c("OBSERVER","DATE","ETIME"))
End.watch=subset(ERSurveyData,subset=EFLAG==5&EXPERIMENT%in%c(1,2),select=c("OBSERVER","DATE","ETIME"))
Start.watch=Start.watch[order(Start.watch$DATE,Start.watch$ETIME),]
End.watch=End.watch[order(End.watch$DATE,End.watch$ETIME),]
RecentExp=Start.watch
RecentExp$Hours=End.watch$ETIME-Start.watch$ETIME
RecentExp=subset(RecentExp,select=c("OBSERVER","DATE","Hours"))
names(RecentExp)[1:2]=c("Observer","Date")
EarlyExp$Observer=as.character(EarlyExp$Observer)
RecentExp$Observer=as.character(RecentExp$Observer)
RecentExp$Date=as.character(RecentExp$Date)
EarlyExp$Date=as.character(EarlyExp$Date)
# Merge data from each set of surveys
ObserverExp=rbind(EarlyExp,RecentExp)
ObserverExp$Observer=factor(ObserverExp$Observer)
ObserverExp$Date=factor(ObserverExp$Date)
# Apply to get sum of effort by observer and date
ObserverExp=t(tapply(ObserverExp$Hours,list(ObserverExp$Observer,ObserverExp$Date),sum))
ObserverExp[is.na(ObserverExp)]=0
# Create cumsum over dates
ObserverExp=apply(ObserverExp,2,cumsum)
# Create a 0 row for first date and then offset dates such that for
# day x the experience is through day x-1
odates=row.names(ObserverExp)
ObserverExp=rbind(rep(0,ncol(ObserverExp)),ObserverExp)
ObserverExp=ObserverExp[-nrow(ObserverExp),]
ObserverExp=as.data.frame(ObserverExp)
ObserverExp$Date=odates
row.names(ObserverExp)=NULL
# Next sum columns for DJR with 35 and CDV with 33 because they had different codes
# amongst the datasets
ObserverExp$DJR=ObserverExp$DJR+ObserverExp[,"35"]
ObserverExp$CDV=ObserverExp$CDV+ObserverExp[,"33"]
return(ObserverExp)
}

