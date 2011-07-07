################################################################################
# Creates dbf files for ERAbund with Recent survey data
################################################################################
CreateERAbundFiles=function(directory="c:/gw",DBDirectory="",years=c(87,92,93,95,97,0,1,6))
{
library(RODBC)
library(foreign)
# Create function to write out the ERAbund dbf file
write.ERAbundFile=function(year,directory=directory)
{
   extract=ERSurveyData[substr(ERSurveyData$SEQUENCE,1,4)==year,2:24]
   extract$DATE=as.Date(as.character(extract$DATE))
   extract$ETIME=as.single(extract$ETIME)
   extract$NTIME=as.single(extract$NTIME)
   extract$STIME=as.single(extract$STIME)
   extract$SRET=as.single(extract$SRET)
   extract$NRET=as.single(extract$NRET)
   extract$C_CPAIR=as.character(extract$C_CPAIR)
   extract=extract[extract$EFLAG!="",]
   write.dbf(extract,paste(directory,"\\/ERSW",substr(year,3,4),substr(year+1,3,4),".dbf",sep=""))
}
# Get recent survey data but first check that it can be found
if(DBDirectory != "")
{
   fdir = file.path(DBDirectory,"GrayWhaleSurveyData.accdb")
} else
{
   fdir = "GrayWhaleSurveyData.accdb"
}
if(!file_test("-f", fdir)) stop(paste("Cannot find file: ",fdir))
con=odbcConnectAccess2007(fdir)
ERSurveyData=sqlFetch(con,"AllRecentData")
# Next check to make sure that all of the ERAbund directories are found and write out the files
for (year in years)
{
  if(directory != "")
  {
   fdir = file.path(directory,paste("ERSW",formatC(year,width=2,flag="0"),formatC(year+1,width=2,flag="0"),sep=""))
  } else
  {
   fdir = file.path(paste("ERSW",formatC(year,width=2,flag="0"),formatC(year+1,width=2,flag="0"),sep=""))
  }
  if(!file_test("-d", fdir)) stop(paste("Cannot find directory: ",fdir))
  fyear=year+1900
  if(fyear<1960)fyear=fyear+100
  write.ERAbundFile(fyear,fdir)
}
# Next do the same for the Secondary subdirectory
for (year in years)
{
  if(directory != "")
  {
   fdir = file.path(directory,paste("Secondary/ERSW",formatC(year,width=2,flag="0"),formatC(year+1,width=2,flag="0"),sep=""))
  } else
  {
   fdir = file.path(paste("Secondary/ERSW",formatC(year,width=2,flag="0"),formatC(year+1,width=2,flag="0"),sep=""))
  }
  if(!file_test("-d", fdir)) stop(paste("Cannot find directory: ",fdir))
  fyear=year+1900
  if(fyear<1960)fyear=fyear+100
  write.ERAbundFile(fyear,fdir)
}
close(con)
return(NULL)
}
# Process ER survey data and create dataframes:
# Primary,PrimaryOff,Secondary,Match,PrimaryEffort,and PrimarySightings
# This function is called by store.ERAbundData which stores the constructed
process.ERdata=function(fdir="c:/gw",maxbeauf=4,maxvis=4)
{
# Get survey data from package directory
  EarlyEffort=NULL
  EarlySightings=NULL
  ERSurveyData=NULL
  data(EarlyEffort,package="ERAnalysis",envir=environment())
  data(EarlySightings,package="ERAnalysis",envir=environment())
  data(ERSurveyData,package="ERAnalysis",envir=environment())
  card2=EarlyEffort
  card3=EarlySightings
  card3$original.watch=card3$Watchperiod
  card3$watch=card3$Assignedwatchperiod
# Read in analysis output files from ERAbund using ReadERAbundFiles
  ERAbundFileList=ReadERAbundFiles(fdir)
  Primary=ERAbundFileList$Primary
  Secondary=ERAbundFileList$Secondary
  Match=ERAbundFileList$Match
  er.migdata=subset(ERAbundFileList$ERNorm,select=c("Start.year","period","begin","end","npods","nwhales","effort","vis","beaufort"))
  er.migdata.sec=subset(ERAbundFileList$ERNormSec,select=c("Start.year","period","begin","end","npods","nwhales","effort","vis","beaufort"))
# Convert some of the fields; distance converted to km
# station is a 5 character string with only last character needed
# watch field created based on crossing time and made into a factor
  Match$distance=Match$distance/1000
  Match$station=factor(substr(as.character(Match$station),5,5))
  Match$original.watch=Match$watch
  Match$watch=1
  Match$watch[Match$t241>=10.5&Match$t241<13.5]=2
  Match$watch[Match$t241>=13.5]=3
  Match$watch=factor(Match$watch)
  Primary$distance=Primary$distance/1000                  # from meters to km
  Primary$original.watch=Primary$watch
  names(Primary)[names(Primary)=="pods.per.hour"]="pphr"
  Primary$etime=NULL
  Primary$station=NULL
  Primary$seen=NULL
  Primary$date=NULL
  Primary$effort.period=paste(Primary$Start.year,"_",formatC(Primary$effort.period,digits=3,flag=0),sep="")
  Primary$watch=1
  Primary$watch[Primary$t241>=10.5]=2
  Primary$watch[Primary$t241>=13.5]=3
  Primary$watch=factor(Primary$watch)
  x=Match[Match$station=="P"&Match$seen==1,]
  x$date=as.Date(paste(x$Start.year,"-12-01",sep=""))+x$day-1
  x$skey=paste(x$date,x$t241)
  Primary$skey=paste(Primary$year,"-",formatC(Primary$month,flag=0,digits=1),
    "-",formatC(Primary$day,flag=0,digits=1)," ",Primary$t241,sep="")
  Primary$only=TRUE
  Primary$only[Primary$skey%in%x$skey]=FALSE
  Primary$skey=NULL
# Create common dataframe for migration curve data
# For early data, create the gwnorm data structure by selecting effort periods
# that meet vis and beaufort criterion
  merged.cards=merge(card2,card3,all.x=TRUE,by.x="key",by.y="key")
# set max beauf and vis values to select data; in data values for Use
# are 4 for both
  badsight=sapply(split(merged.cards,merged.cards$key), function(x) any(x$Visibility>maxvis | x$Beaufort>maxbeauf))
  badsight[is.na(badsight)]=FALSE
  bb=subset(card2,select=c("Block1beaufort","Block1visibility","Block2beaufort","Block2visibility","Block3beaufort","Block3visibility"))
  badbf=apply(bb[,c(1,3,5)],1,function(x)any(x>maxbeauf))
  badvis=apply(bb[,c(2,4,6)],1,function(x)any(x>maxvis))
  badbf[is.na(badbf) | is.nan(badbf)]=FALSE
  badvis[is.na(badvis)| is.nan(badvis)]=FALSE
  card2$Use="Yes"
  card2$Use[badvis | badbf | badsight]="No"
# create an average value for vis and beaufort using the sightings and effort values for the older data
  avgviseff=apply(subset(card2,select=c("Block1visibility","Block2visibility","Block3visibility")),1,mean,na.rm=TRUE)
  avgvis=with(merged.cards, tapply(Visibility,key,mean,na.rm=TRUE))
  avgvis[is.na(avgvis)]=avgviseff[is.na(avgvis)]
  avgbfeff=apply(subset(card2,select=c("Block1beaufort","Block2beaufort","Block3beaufort")),1,mean,na.rm=TRUE)
  avgbf=with(merged.cards, tapply(Beaufort,key,mean,na.rm=TRUE))
  avgbf[is.na(avgbf)]=avgbfeff[is.na(avgbf)]
  card2$vis=as.vector(avgvis)
  card2$beaufort=as.vector(avgbf)
# merge effort and sightings again with new Use values
  merged.cards=merge(card2,card3,all.x=TRUE,by.x="key",by.y="key")
# construct migration effort dataframe for early years
  Use=NULL
  gwnorm=subset(card2,subset=Use=="Yes",select=c("key","Start.year","Year","Month","Day","Begin.date.time","End.date.time","vis","beaufort","Observer"))
# construct a dataframe of pod and whale counts from observations
  npods=with(merged.cards[merged.cards$Use=="Yes" &merged.cards$Included,], tapply(Podsize,key,length))
  nwhales=with(merged.cards[merged.cards$Use=="Yes" &merged.cards$Included,], tapply(Podsize,key,sum))
  df=data.frame(key=names(npods),npods=as.vector(npods),nwhales=as.vector(nwhales))
# merge the counts with the effort data
  gwnorm=merge(gwnorm,df,all.x=TRUE,by="key")
  gwnorm$npods[is.na(gwnorm$npods)]=0
  gwnorm$nwhales[is.na(gwnorm$nwhales)]=0
# compute length of effort period in decimal days and compute begin and
# end values for period in decimal days from 1 Dec
  gwnorm$begin=gwnorm$Begin.date.time-strptime(paste(gwnorm$Start.year,"-12-01 00:00",sep=""),format="%Y-%m-%d %H:%M")
  gwnorm$end=gwnorm$End.date.time-strptime(paste(gwnorm$Start.year,"-12-01 00:00",sep=""),format="%Y-%m-%d %H:%M")
  gwnorm$effort=as.vector(gwnorm$end)-as.vector(gwnorm$begin)
  gwnorm=subset(gwnorm,select=c("Start.year","key","begin","end","npods","nwhales","effort","vis","beaufort","Observer"))
# Create single primary sightings file with data from all years
  PrimarySightings=merge(subset(card3,select=-Observer),gwnorm,by.x="key",by.y="key")
  PrimarySightings$only=TRUE
  PrimarySightings=PrimarySightings[PrimarySightings$Included,]
  PrimarySightings$distance=PrimarySightings$Distanceoffshore*1.852 # from nm to km
  PrimarySightings$distance[is.na(PrimarySightings$distance)]=mean(PrimarySightings$distance[!is.na(PrimarySightings$distance)])
  PrimarySightings$pphr=PrimarySightings$npods/(PrimarySightings$effort*24)
  names(PrimarySightings)[10]="podsize"
  names(PrimarySightings)[15]="vis"
  names(PrimarySightings)[31]="Start.year"
  T241c=as.character(PrimarySightings$Timeabeam)
  PrimarySightings$t241=as.numeric(substr(T241c,1,2))+as.numeric(substr(T241c,4,5))/60+as.numeric(substr(T241c,7,8))/3600
  PrimarySightings=subset(PrimarySightings,select=c("Day","Month","Year","watch","t241","distance",
             "podsize","Observer","vis","Beaufort","Winddirection","key","pphr","Start.year","original.watch","only"))
  names(PrimarySightings)=names(Primary)
  PrimarySightings$pphr=as.numeric(PrimarySightings$pphr)
  PrimarySightings=rbind(PrimarySightings,Primary)
  names(PrimarySightings)[names(PrimarySightings)=="effort.period"]="key"
  names(Primary)[names(Primary)=="effort.period"]="key"
  PrimarySightings$watch=factor(PrimarySightings$watch)
# merge newer effort with older effort in gwnorm
  names(er.migdata)[names(er.migdata)=="period"]="key"
  er.migdata$key=paste(er.migdata$Start.year,"_",formatC(er.migdata$key,digits=3,flag=0),sep="")
  gwnorm$key=as.character(gwnorm$key)
#  er.migdata=rbind(subset(gwnorm,select=-Observer),er.migdata)
  er.migdata$Observer=0
  er.migdata=rbind(gwnorm,er.migdata)
  er.migdata$time=as.vector((er.migdata$end+er.migdata$begin)/2)
  er.migdata$begin=as.numeric(er.migdata$begin)
  er.migdata$end=as.numeric(er.migdata$end)
  er.migdata$key=factor(er.migdata$key)
# define watch 1-3 to select data based on
# the presence of any beauf>maxbeauf or vis>maxvis as in earlier years
  er.migdata$watch=1
  er.migdata$watch[er.migdata$Start.year<1987][(er.migdata$time[er.migdata$Start.year<1987]-floor(er.migdata$time[er.migdata$Start.year<1987]))>=(12/24)]=2
  er.migdata$watch[er.migdata$Start.year>=1987][(er.migdata$time[er.migdata$Start.year>=1987]-floor(er.migdata$time[er.migdata$Start.year>=1987]))>=(10.5/24)]=2
  er.migdata$watch[er.migdata$Start.year>=1987][(er.migdata$time[er.migdata$Start.year>=1987]-floor(er.migdata$time[er.migdata$Start.year>=1987]))>=(13.5/24)]=3
  er.migdata$watch.key=paste(er.migdata$Start.year,"_",floor(er.migdata$begin),er.migdata$watch,sep="")
  notUse=sapply(split(er.migdata,er.migdata$watch.key),function(x) any(x$vis>maxvis) | any(x$beaufort>maxbeauf))
  notUse[is.na(notUse)]=TRUE
  er.migdata=merge(er.migdata,data.frame(watch.key=names(notUse),Use=!notUse))
#
# Restrict all PrimarySightings & PrimaryEffort to be within maxbeauf and maxvis
  PrimaryEffort=er.migdata[(er.migdata$vis<=maxvis| is.na(er.migdata$vis))&(is.na(er.migdata$beaufort) | er.migdata$beaufort<=maxbeauf),]
  PrimarySightings=PrimarySightings[(is.na(PrimarySightings$vis)|PrimarySightings$vis<=maxvis)&(is.na(PrimarySightings$beaufort)|PrimarySightings$beaufort<=maxbeauf),]
# Create observer experience dataframe and add the experience field in hours to Match and
# PrimarySightings
  ObserverExp=CreateObserverExperience()
  data(Observer,package="ERAnalysis",envir=environment())
  PrimarySightings$Date=paste(PrimarySightings$year,"-",formatC(PrimarySightings$month,width=2,flag=0),"-",formatC(PrimarySightings$day,width=2,flag=0),sep="")
  PrimarySightings$observer=substr(PrimarySightings$observer,nchar(PrimarySightings$observer)-2,nchar(PrimarySightings$observer))
  Observer$key=as.character(Observer$Initials)
  Observer$key[is.na(Observer$Initials)]=Observer$Observer[is.na(Observer$Initials)]
  PrimarySightings=merge(PrimarySightings,ObserverExp,by="Date",all.x=TRUE)
  PrimarySightings$hours=PrimarySightings[cbind(1:nrow(PrimarySightings),match(PrimarySightings$observer,names(PrimarySightings)[18:76])+17)]
  PrimarySightings=PrimarySightings[,c(1:17,77)]
  PrimarySightings$seq=1:nrow(PrimarySightings)
  PrimarySightings=merge(PrimarySightings, subset(Observer, select = c("key", "Sex")), by.x ="observer",by.y = "key", all.x = TRUE)
  PrimarySightings=PrimarySightings[order(PrimarySightings$seq),]
  PrimarySightings$seq=NULL
  PrimarySightings$hours=as.numeric(PrimarySightings$hours)
  PrimaryEffort$Date=substr(as.character(strptime(paste(PrimaryEffort$Start.year,"-12-01",sep=""),format="%Y-%m-%d")+3600*24*(floor(PrimaryEffort$begin)-1)),1,10)
  PrimarySightings$observer[PrimarySightings$observer=="33"]="CDV"
  PrimarySightings$observer[PrimarySightings$observer=="35"]="DJR"
  PrimarySightings$Observer=factor(PrimarySightings$observer)
  PrimarySightings$observer=NULL
  Match$Date=as.character(strptime(paste(Match$Start.year,"-12-01",sep=""),format="%Y-%m-%d")+3600*24*(Match$day-1))
  Match=merge(Match,ObserverExp,by="Date",all.x=TRUE)
  Match$observer=as.character(Match$observer)
  Match$observer=substr(Match$observer,nchar(Match$observer)-2,nchar(Match$observer))
  Match$hours=Match[cbind(1:nrow(Match),match(Match$observer,names(Match)[18:76])+17)]
  Match=Match[,c(1:17,77)]
  Match$hours=as.numeric(Match$hours)
  Match$seq=1:nrow(Match)
  Match=merge(Match, subset(Observer, select = c("key", "Sex")), by.x ="observer",by.y = "key", all.x = TRUE)
  Match=Match[order(Match$seq),]
  Match$seq=NULL
  Match$observer[Match$observer=="33"]="CDV"
  Match$observer[Match$observer=="35"]="DJR"
  Match$Observer=factor(Match$observer)
  Match$observer=NULL
# Create Secondary effort file
  names(er.migdata.sec)[names(er.migdata.sec)=="period"]="key"
  er.migdata.sec$key=paste(er.migdata.sec$Start.year,"_",formatC(er.migdata.sec$key,digits=3,flag=0),sep="")
  er.migdata.sec$time=as.vector((er.migdata.sec$end+er.migdata.sec$begin)/2)
  er.migdata.sec$begin=as.numeric(er.migdata.sec$begin)
  er.migdata.sec$end=as.numeric(er.migdata.sec$end)
  er.migdata.sec$key=factor(er.migdata.sec$key)
# define watch 1-3 to select data based on
# the presence of any beauf>maxbeauf or vis>maxvis as in earlier years
  er.migdata.sec$watch=1
  er.migdata.sec$watch[er.migdata.sec$Start.year<1987][(er.migdata.sec$time[er.migdata.sec$Start.year<1987]-floor(er.migdata.sec$time[er.migdata.sec$Start.year<1987]))>=(12/24)]=2
  er.migdata.sec$watch[er.migdata.sec$Start.year>=1987][(er.migdata.sec$time[er.migdata.sec$Start.year>=1987]-floor(er.migdata.sec$time[er.migdata.sec$Start.year>=1987]))>=(10.5/24)]=2
  er.migdata.sec$watch[er.migdata.sec$Start.year>=1987][(er.migdata.sec$time[er.migdata.sec$Start.year>=1987]-floor(er.migdata.sec$time[er.migdata.sec$Start.year>=1987]))>=(13.5/24)]=3
  er.migdata.sec$watch.key=paste(er.migdata.sec$Start.year,"_",floor(er.migdata.sec$begin),er.migdata.sec$watch,sep="")
  notUse=sapply(split(er.migdata.sec,er.migdata.sec$watch.key),function(x) any(x$vis>maxvis) | any(x$beaufort>maxbeauf))
  er.migdata.sec=merge(er.migdata.sec,data.frame(watch.key=names(notUse),Use=!notUse))
#
# Restrict all SecondaryEffort to be within maxbeauf and maxvis
  SecondaryEffort=er.migdata.sec[(er.migdata.sec$vis<=maxvis| is.na(er.migdata.sec$vis))&(is.na(er.migdata.sec$beaufort) | er.migdata.sec$beaufort<=maxbeauf),]
  SecondaryEffort$Date=substr(as.character(strptime(paste(SecondaryEffort$Start.year,"-12-01",sep=""),format="%Y-%m-%d")+3600*24*(floor(SecondaryEffort$begin)-1)),1,10)
  convert.decimal.time=function(x)
  {
  xx=(x-floor(x))*60
  zz=trunc((xx-floor(xx))*60+.5)
  return(paste(formatC(floor(x),width=2,flag=0),formatC(floor(xx),width=2,flag=0),formatC(floor(zz),width=2,flag=0),sep=":"))
  }
  data(ERSurveyData,package="ERAnalysis",envir=environment())
  nn=names(PrimaryEffort)
# 1987 surveys and beyond -- adding Observer to Primary & Secondary Effort data
# Primary Effort EFLAG=1,2 records
  EFLAG=NULL
  EXPERIMENT=NULL
  LOCATION=NULL
  Start.year=NULL
  xx=subset(ERSurveyData,subset=EFLAG<=2&(EXPERIMENT==1 | (EXPERIMENT==2 & ((LOCATION=="S" & !Start.year %in% 2000:2001)|(Start.year %in% 2000:2001 & LOCATION=="N")) )))
  xx=xx[xx$VISCODE<=maxvis&xx$WINDFORCE<=maxbeauf,]
  xx$Start.year=as.numeric(xx$Start.year)
  xx$key=paste(xx$DATE,"_",substr(convert.decimal.time(xx$ETIME),1,2),sep="")
  PrimaryEffort$seq=1:nrow(PrimaryEffort)
  PrimaryEffort$wkey=paste(PrimaryEffort$Date,substr(convert.decimal.time(24*(PrimaryEffort$begin-floor(PrimaryEffort$begin)+.5/3600)),1,2),sep="_")
  zz=merge(PrimaryEffort,subset(xx,select=c("key","OBSERVER","ETIME")),by.x="wkey",by.y="key",all.x=TRUE)
  zz$tt=24*(zz$begin-floor(zz$begin))
  zz=zz[abs(zz$tt-zz$ETIME)<=1.5/3600,]
  zz=zz[!is.na(zz$wkey),]
  zz=zz[order(zz$seq),]
  PrimaryEffort$Observer=as.character(PrimaryEffort$Observer)
  PrimaryEffort$Observer[PrimaryEffort$Observer==0]=as.character(zz$OBSERVER)
  PrimaryEffort$Observer=factor(PrimaryEffort$Observer)
  PrimaryEffort=subset(PrimaryEffort,select=nn)
# Secondary Effort EFLAG=1,2 records
  xx=subset(ERSurveyData,subset=EFLAG<=2& (EXPERIMENT==2 & ((LOCATION=="S" & Start.year %in% 2000:2001)|(!Start.year %in% 2000:2001 & LOCATION=="N")) ))
  xx=xx[xx$VISCODE<=maxvis&xx$WINDFORCE<=maxbeauf,]
  xx$Start.year=as.numeric(xx$Start.year)
  xx$key=paste(xx$DATE,"_",substr(convert.decimal.time(xx$ETIME),1,2),sep="")
  SecondaryEffort$seq=1:nrow(SecondaryEffort)
  SecondaryEffort$wkey=paste(SecondaryEffort$Date,substr(convert.decimal.time(24*(SecondaryEffort$begin-floor(SecondaryEffort$begin)+.5/3600)),1,2),sep="_")
  zz=merge(SecondaryEffort,subset(xx,select=c("key","OBSERVER","ETIME")),by.x="wkey",by.y="key",all.x=TRUE)
  zz$tt=24*(zz$begin-floor(zz$begin))
  zz=zz[abs(zz$tt-zz$ETIME)<=1.5/3600,]
  zz$Observer=zz$OBSERVER
  SecondaryEffort=zz
  SecondaryEffort=subset(SecondaryEffort,select=nn)
  Primary$Observer=factor(Primary$observer)
  Primary$observer=NULL
  return(list(Primary=Primary,PrimaryOff=ERAbundFileList$PrimaryOff,Secondary=Secondary,Match=Match,PrimarySightings=PrimarySightings,
               PrimaryEffort=PrimaryEffort,SecondaryEffort=SecondaryEffort))
}
################################################################################
# Reads in output data files created by ERAbund with recent survey data
# The function returns a list of 4 dataframes:
#  1) Primary   - on-effort sightings from primary platform
#  2) PrimaryOff  - off-effort sightings from primary platform
#  3) Secondary - on and off-effort sightings from secondary platform
#  4) Match     - matches for double-observer analysis
#  5) ERNorm    - data structure for GWNorm program
################################################################################
ReadERAbundFiles=function(directory="c:/gw",years=c(87,92,93,95,97,0,1,6))
{
Primary=NULL
PrimaryOff=NULL
Secondary=NULL
Match=NULL
ERNorm=NULL
ERNormSec=NULL
for (year in years)
{
  yearc=paste(formatC(year,width=2,flag="0"),formatC(year+1,width=2,flag="0"),sep="")
# Read in primary observation file
  if(directory != "")
  {
   fdir = file.path(directory,paste("ERSW",yearc,sep=""),paste("erprim",yearc,".dat",sep=""))
  } else
  {
   fdir = file.path(paste("ERSW",yearc,sep=""),paste("erprim",yearc,".dat",sep=""))
  }
  if(!file_test("-f", fdir)) stop(paste("Cannot find directory: ",fdir))
  Primaryx=read.fwf(fdir,widths=c(5,5,5,9,4,11,10,4,5,4,5,5,8,7,5,5),skip=1)
  names(Primaryx)=c("day","month","year","etime","watch","t241","distance","podsize",
       "observer","vis","beaufort","wind.direction","effort.period","pods.per.hour","station","seen")
  Primaryx$Start.year=year+1900
  Primaryx$Start.year[Primaryx$Start.year<1987]=Primaryx$Start.year[Primaryx$Start.year<1987]+100
  Primary=rbind(Primary,Primaryx)
# Read in primary observation 'not used' file
  if(directory != "")
  {
   fdir = file.path(directory,paste("ERSW",yearc,sep=""),paste("erprimUAE",yearc,".dat",sep=""))
  } else
  {
   fdir = file.path(paste("ERSW",yearc,sep=""),paste("erprimUAE",yearc,".dat",sep=""))
  }
  if(!file_test("-f", fdir)) stop(paste("Cannot find directory: ",fdir))
  Primaryx=read.fwf(fdir,widths=c(5,5,5,9,4,11,10,4,5,4,5,5,8,7,5,5),skip=1)
  names(Primaryx)=c("day","month","year","etime","watch","t241","distance","podsize",
       "observer","vis","beaufort","wind.direction","effort.period","pods.per.hour","station","seen")
  Primaryx$Start.year=year+1900
  Primaryx$Start.year[Primaryx$Start.year<1987]=Primaryx$Start.year[Primaryx$Start.year<1987]+100
  PrimaryOff=rbind(PrimaryOff,Primaryx)
# Read in secondary observation file
  if(directory != "")
  {
   fdir = file.path(directory,paste("ERSW",yearc,sep=""),paste("ersec",yearc,".dat",sep=""))
  } else
  {
   fdir = file.path(paste("ERSW",yearc,sep=""),paste("ersec",yearc,".dat",sep=""))
  }
  if(!file_test("-f", fdir)) stop(paste("Cannot find directory: ",fdir))
  Secondaryx=read.fwf(fdir,widths=c(5,5,5,9,11,12,4,4,3,5,10))
  names(Secondaryx)=c("day","month","year","etime","t241","distance","podsize",
       "vis","beaufort","wind.direction","off")
  Secondaryx$Start.year=year+1900
  Secondaryx$Start.year[Secondaryx$Start.year<1987]=Secondaryx$Start.year[Secondaryx$Start.year<1987]+100
  Secondary=rbind(Secondary,Secondaryx)
# Read in match observation file
  if(directory != "")
  {
   fdir = file.path(directory,paste("ERSW",yearc,sep=""),paste("ermatchsp",yearc,".dat",sep=""))
  } else
  {
   fdir = file.path(paste("ERSW",yearc,sep=""),paste("ermatchsp",yearc,".dat",sep=""))
  }
  if(!file_test("-f", fdir)) stop(paste("Cannot find directory: ",fdir))
  Matchx=read.fwf(fdir,widths=c(5,5,5,4,9,10,6,4,4,5,5,8,9,6),skip=1)
  names(Matchx)=c("seen","station","day","watch","t241","distance","podsize","observer",
       "vis","beaufort","wind.direction","pphr","mscore","mcode")
  Matchx$Start.year=year+1900
  Matchx$Start.year[Matchx$Start.year<1987]=Matchx$Start.year[Matchx$Start.year<1987]+100
  Match=rbind(Match,Matchx)
# Read in ernorm file
  if(directory != "")
  {
   fdir = file.path(directory,paste("ERSW",yearc,sep=""),paste("ernorm",yearc,".dat",sep=""))
  } else
  {
   fdir = file.path(paste("ERSW",yearc,sep=""),paste("ernorm",yearc,".dat",sep=""))
  }
  if(!file_test("-f", fdir)) stop(paste("Cannot find directory: ",fdir))
  normx=read.fwf(fdir,widths=c(7,12,11,4,4,rep(3,13),4))
  names(normx)=c("period","begin","end","npods","nwhales",paste("pod",c(1:10,"11+"),sep=""),"vis","beaufort","wind.direction")
  normx$Start.year=year+1900
  normx$Start.year[normx$Start.year<1987]=normx$Start.year[normx$Start.year<1987]+100
  ERNorm=rbind(ERNorm,normx)
  # Read in secondary ernorm file
  if(directory != "")
  {
   fdir = file.path(directory,paste("Secondary/ERSW",yearc,sep=""),paste("ernorm",yearc,".dat",sep=""))
  } else
  {
   fdir = file.path(paste("Secondary/ERSW",yearc,sep=""),paste("ernorm",yearc,".dat",sep=""))
  }
  if(!file_test("-f", fdir)) stop(paste("Cannot find directory: ",fdir))
  normx=read.fwf(fdir,widths=c(7,12,11,4,4,rep(3,13),4))
  names(normx)=c("period","begin","end","npods","nwhales",paste("pod",c(1:10,"11+"),sep=""),"vis","beaufort","wind.direction")
  normx$Start.year=year+1900
  normx$Start.year[normx$Start.year<1987]=normx$Start.year[normx$Start.year<1987]+100
  ERNormSec=rbind(ERNormSec,normx)
}
Primary$wind.direction=factor(Primary$wind.direction)
Primary$year=Primary$year+1900
Primary$year[Primary$year<1987]=Primary$year[Primary$year<1987]+100
Primary$date=as.Date(paste(Primary$year,Primary$month,Primary$day,sep="-"))
Secondary$wind.direction=factor(Secondary$wind.direction)
Secondary$year=Secondary$year+1900
Secondary$year[Secondary$year<1987]=Secondary$year[Secondary$year<1987]+100
Secondary$date=as.Date(paste(Secondary$year,Secondary$month,Secondary$day,sep="-"))
Match$observer=factor(Match$observer)
Match$wind.direction=factor(Match$wind.direction)
ERNorm$effort=as.vector(ERNorm$end-ERNorm$begin)
ERNormSec$effort=as.vector(ERNormSec$end-ERNormSec$begin)
return(list(Primary=Primary,PrimaryOff=PrimaryOff,Secondary=Secondary,Match=Match,ERNorm=ERNorm,ERNormSec=ERNormSec))
}
store.ERdata=function(DBDirectory="",package.dir="C:/Users/Jeff Laake/Desktop/MyDocuments/R Development/ERAnalysis/ERAnalysis/data",nmax=20)
{
# Retrieves data from Gray Whale Access Database and creates dataframes
# EarlyEffort, EarlySightings for surveys from 1967-1985 and ERSurveyData for
# surveys from 1987-2006.  Some computed fields are added to the dataframes and then
# they are saved to the ERAnalysis package data directory.
  library(RODBC)
  # Set DBDirectory to point to Access Database. "" means in local directory
  DBDirectory=""
  if(DBDirectory != "")
  {
     fdir = file.path(DBDirectory,"GrayWhaleSurveyData.accdb")
  } else
  {
     fdir = "GrayWhaleSurveyData.accdb"
  }
  if(!file_test("-f", fdir)) stop(paste("Cannot find file: ",fdir))
  con=odbcConnectAccess2007(fdir)
# get effort and sightings tables and save them in package directory
  EarlyEffort=sqlFetch(con,"EffortPre1987")
  EarlySightings=sqlFetch(con,"ERSightingsPre1987")
  ERSurveyData=sqlFetch(con,"AllRecentData")
  EarlyEffort$Begin.date.time=as.POSIXct(with(EarlyEffort,strptime(paste(Year,"-",Month,"-",Day," ",formatC(Begintime,digits=3,flag="0"),sep=""),
                                  format="%Y-%m-%d %H%M")))
  EarlyEffort$End.date.time=as.POSIXct(with(EarlyEffort,strptime(paste(Year,"-",Month,"-",Day," ",formatC(Endtime,digits=3,flag="0"),sep=""),
                                  format="%Y-%m-%d %H%M")))
  EarlyEffort$Start.year=EarlyEffort$Year
  EarlyEffort$Start.year[EarlyEffort$Month<4]=EarlyEffort$Start.year[EarlyEffort$Month<4]-1
  EarlySightings$Start.year=EarlySightings$Year
  EarlySightings$Start.year[EarlySightings$Month<4]=EarlySightings$Start.year[EarlySightings$Month<4]-1
  ERSurveyData$Start.year=substr(ERSurveyData$SEQUENCE,1,4)
  ERSurveyData$year=substr(ERSurveyData$DATE,1,4)
  ERSurveyData$month=substr(ERSurveyData$DATE,6,7)
  ERSurveyData$day=substr(ERSurveyData$DATE,9,10)
  save(EarlyEffort,file=file.path(package.dir,"EarlyEffort.rda"))
  save(EarlySightings,file=file.path(package.dir,"EarlySightings.rda"))
  save(ERSurveyData,file=file.path(package.dir,"ERSurveyData.rda"))
  Observer=sqlFetch(con,"ObserverCodes")
  AddObs=subset(Observer,subset=Observer%in%c(33,35))
  AddObs$Initials=NA
  Observer=rbind(Observer,AddObs)
  save(Observer,file=file.path(package.dir,"Observer.rda"))
  close(con)
  # Read in calibration data from text file and create and store dataframes
  PodsizeCalibrationData=read.delim(file.path(package.dir,"PScaliball.txt"),header=TRUE)
  PodsizeCalibrationData$key=paste(PodsizeCalibrationData$Type,PodsizeCalibrationData$Year,PodsizeCalibrationData$ID,sep="_")
  PodsizeCalibrationTable=as.data.frame(table(PodsizeCalibrationData$key,factor(PodsizeCalibrationData$Estimate,levels=1:nmax))[,1:nmax])
  PodsizeCalibrationTable$key=row.names(PodsizeCalibrationTable)
  row.names(PodsizeCalibrationTable)=NULL
  PodsizeCalibrationTable=merge(unique(data.frame(key=PodsizeCalibrationData$key,
        True=PodsizeCalibrationData$True)),PodsizeCalibrationTable,by="key")
  save(PodsizeCalibrationTable,file=file.path(package.dir,"PodsizeCalibrationTable.rda"))
  save(PodsizeCalibrationData,file=file.path(package.dir,"PodsizeCalibrationData.rda"))
}
store.ERAbundData=function(fdir="c:/gw",package.dir="C:/Users/Jeff Laake/Desktop/MyDocuments/R Development/ERAnalysis/ERAnalysis/data",
                                 maxbeauf=4,maxvis=4)
{
   df.list=process.ERdata(fdir=fdir,maxbeauf=maxbeauf,maxvis=maxvis)
   SecondarySightings=df.list$Secondary
   Primary=df.list$Primary
   PrimarySightings=df.list$PrimarySightings
   PrimaryEffort=df.list$PrimaryEffort
   SecondaryEffort=df.list$SecondaryEffort
   Match=df.list$Match
   PrimaryOff=df.list$PrimaryOff
   save(Primary,file=file.path(package.dir,"Primary.rda"))
   save(PrimarySightings,file=file.path(package.dir,"PrimarySightings.rda"))
   save(PrimaryEffort,file=file.path(package.dir,"PrimaryEffort.rda"))
   save(SecondaryEffort,file=file.path(package.dir,"SecondaryEffort.rda"))
   save(Match,file=file.path(package.dir,"Match.rda"))
   save(SecondarySightings,file=file.path(package.dir,"SecondarySightings.rda"))
   save(PrimaryOff,file=file.path(package.dir,"PrimaryOff.rda"))
   df.list=process.ERdata(fdir=fdir,maxbeauf=10,maxvis=10)
   AllPrimaryEffort=df.list$PrimaryEffort
   AllSecondaryEffort=df.list$SecondaryEffort
   save(AllPrimaryEffort,file=file.path(package.dir,"AllPrimaryEffort.rda"))
   save(AllSecondaryEffort,file=file.path(package.dir,"AllSecondaryEffort.rda"))
}



