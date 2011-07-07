##########################################################################
#  Function linklist links primary and secondary observer records in Match.list
#  creates revised dataframe with t241 and distance averaged and group size summed
#  three possible types of spatial criterion are included:
# "box" criterion is maximum of distances
# "diamond" criterion is sum of distances
# "ellipse" criterion is the euclidean distance default is ellipse.
#  Three specific criteria based on prior analysis are included:
#  "ERAbund" Uses the diamond spatial criterion with user specified parameters,
#         default is the criteria used by Breiwick et all 2006
# 0.04   !Time weight factor for linking
# 1.75   !Distance weight factor for linking
# 0.00   !Pod size weight factor for linking
#-0.05   !Maximum sum score for linking (merging pods)
## Rugh et al. 1993
# Rugh et al. 1990
#  RCH 6/12/2009
##################################################################
linklist <- function(mspdf,lcut=-0.05,twt=.18,dwt=3.95,crittype="ellipse")
{
# Arguments:
#
# mspdf  - dataframe of sightings
# t2dist - multiplier to convert time to distancein km/hr
# crittype - criterion for linking
# lcut   - linking parameter cutoff
# twt  -time weight for crossing time difference in minutes
# dwt   distance weight for ratio of difference to longer distance of shore at 241
#
mdf <- as.data.frame(mspdf)
Pdf <- mdf[(mdf$station=="P" & mdf$t241>1 ),]
Sdf <- mdf[(mdf$station=="S" & mdf$t241>1 ),]
lenp <- nrow(Pdf)
lens <- nrow(Sdf)
if (lenp>=2 & lcut>= 0) {
    Pdf <- linkdf(Pdf,lcut,twt,dwt,crittype)
    lenp <- nrow(Pdf)
    }
if (lens>=2 & lcut>= 0) {
    Sdf <- linkdf(Sdf,lcut,twt,dwt,crittype)
    lens <- nrow(Sdf)
    }
linklist<-rbind(Pdf,Sdf)
return(linklist)
}

##########################################################################
#  Function linkdf identifies potential links within observer record and
#  creates revised dataframe with t241 and dist averaged and group size summed
#  three possible types of spatial criterion are included:
# "box" criterion is maximum of distances
# "diamond" criterion is sum of distances
# "ellipse" criterion is the euclidean distance default is ellipse.
#  Three specific criteria based on prior analysis are included:
#  "ERAbund" Uses the diamond spatial criterion with user specified parameters,
#         default is the criteria used by Breiwick et all 2006
# 0.04   !Time weight factor for linking
# 1.75   !Distance weight factor for linking
# 0.00   !Pod size weight factor for linking
#-0.05   !Maximum sum score for linking (merging pods)
# Rugh et al. 1993
# Rugh et al. 1990
#  RCH 6/12/2009
####################################################################
linkdf <- function(mspdf,lcut=-0.05,twt=.18,dwt=3.95,crittype="ellipse")
{
# Arguments:
#
# mspdf  - dataframe of sightings
# t2dist - multiplier to convert time to distancein km/hr
# crittype - criterion for linking
# lcut   - linking parameter cutoff
# twt  -time weight for crossing time difference in minutes
# dwt   distance weight for ratio of difference to longer distance of shore at 241
#
mdf <- as.data.frame(mspdf)
lenm <- nrow(mdf)
if (lenm<=1) return(mdf)
tmat <- twt*60*abs(outer(mdf$t241, mdf$t241, "-"))
 fd <- function(x, y) (x-y)/max(x,y)
 dmat <- dwt*abs(outer(mdf$distance, mdf$distance, fd))
if (crittype == "box")   cmat<- pmax(tmat,dmat)
else
  if (crittype == "diamond")   cmat<- tmat + dmat
else
  if (crittype == "ERAbund")     cmat<- tmat + dmat
else   cmat<- sqrt(tmat^2 + dmat^2)
#generate link matrix
 c1mat <- 1*( cmat < lcut & cmat > 0)
#generate conversion matrix
 comat<-c1mat + diag(lenm)
#count entries byrow
crvec<-  rowSums(comat)
ccvec<-  colSums(comat)
#weight columns
conmat<-  comat * sqrt(outer(1/crvec,1/ccvec))
#sum weights by row
 smvec <-  rowSums(conmat)
#check sum if all rows are prperly weighted then finished
 while(length(smvec[smvec != 1])>0)
 {
#find improperly weighted rows
   snvec<- (smvec!=1)
#mask matrix to get subset of rows to be fixed
   fmat<-outer(snvec,snvec)  *cmat * c1mat
 # find and remove maximum link distance
   c1mat<-c1mat * (fmat!= max(fmat))
# create new conversion matrix
   comat<-c1mat + diag(lenm)
crvec<-  rowSums(comat)
ccvec<-  colSums(comat)
conmat<-  comat * sqrt(outer(1/crvec,1/ccvec))
  #check sum
     smvec <-  rowSums(conmat)
 }
#find rows with elements in lower triangle indicating duplicates
convec<- rowSums(lower.tri(conmat)*conmat)
# generate linked db
linkdf<-mdf
#t241 and dist averaged
linkdf$t241<- as.vector(conmat%*%(mdf$t241 ))
linkdf$distance<- as.vector(conmat%*%(mdf$distance ))
#ps summed
sumps<- (1*(conmat>0))
linkdf$podsize<- as.vector(sumps%*%(mdf$podsize))
if(!is.null(mdf$corrected.podsize))linkdf$corrected.podsize<- as.vector(sumps%*%(mdf$corrected.podsize))
#remove duplicates
linkdf<- linkdf[ convec==0,]
return(linkdf)
}

##########################################################################
#  Function matchdf identifies potential links within a linked observer record and
#  creates a revised dataframe with matches ordered by primary then secondary
#  column 1 is seen column 2 is 1 for primary and 2 for secondary
#  three possible types of spatial criterion are included:
# "box" criterion is maximum of distances
# "diamond" criterion is sum of distances
# "ellipse" criterion is the euclidean distance default is ellipse.
#  Three specific criteria based on prior analysis are included:
#  "ERAbund" Uses the diamond spatial criterion with user specified parameters,
#         default is the criteria used by Breiwick et all 200??
#0.04         !Time weight factor for matching     time difference in minutes
#1.75         !Distance weight factor for matching
#0.05         !Pod size weight factor for matching
# 1.0         !Maximum sum score for matching
# Rugh et al. 1993
# Rugh et al. 1990
#  RCH 7/21/2009
####################################################################

 matchdf<- function(mspdf,mcut=1,twt=.18,dwt=3.95,pwt=.05,crittype="ellipse")
{
# Arguments:
#
# mspdf  - dataframe of sightings
# t2dist - multiplier to convert time to distancein km/hr
# crittype - criterion for matching
# mcut   - matching parameter cutoff
# twt  -time weight for crossing time difference in minutes
# dwt   distance weight for ratio of difference to longer distance of shore at 241
# pwt   weight for difference in podsize
#
mdf <- as.data.frame(mspdf)
Pdf <- mdf[(mdf$station=="P"),]
Sdf <- mdf[(mdf$station=="S"),]
lenp <- nrow(Pdf)
lens <- nrow(Sdf)
if((lenp*lens)>0)  {
tmat <- twt*60*abs(outer(Pdf$t241, Sdf$t241, "-"))
 fd <- function(x, y) (x-y)/max(x,y)
 dmat <- dwt*abs(outer(Pdf$distance, Sdf$distance, fd))
 pmat <- pwt*abs(outer(Pdf$podsize, Sdf$podsize, "-"))
if (crittype == "box")     cmat<- pmax(tmat,dmat) + pmat
else
  if (crittype == "diamond")     cmat<- tmat + dmat + pmat
else
  if (crittype == "ERAbund")     cmat<- tmat + dmat + pmat
  else     cmat<- sqrt(tmat^2 + dmat^2) + pmat
# sort by best and eliminate sequentialy from best to worst
omat<- cbind(c(cmat),rep(1:lenp,lens),rep(1:lens,each=lenp))
      ii<-order(omat[,1],omat[,2])
     omat1<-rbind(omat[ii,],9999999,9999999)
     matched<-rbind(c(0,0,0),c(0,0,0))
     while(omat1[1,1] <= mcut){
     om2<-  omat1[1,2]
     om3<-  omat1[1,3]
     matched<-rbind(matched,c(1,omat1[1,2:3] ))
     omat1<-omat1[omat1[,2]!= om2,]
     omat1<-omat1[omat1[,3]!= om3,]
       } #end while
# nomatch vectors
crvec<- (1:lenp %in% matched[,2])                       
ccvec<- (1:lens %in% matched[,3])                       
xr<-cbind(crvec,2,1:lenp,0)
xc<-cbind(ccvec,3,0,1:lens)
#build matchdf
matches<-rbind(matched[ matched[,1]==1,],xr[xr[,1]==0,2:4],xc[xc[,1]==0,2:4])
} # end of if((lenp*lens)>0)
# generate matches set if no sightings by one observer
else {
if (lenp >0) matches<-cbind(2,1:lenp,0)
else matches<-cbind(3,0,1:lens)
}   #end else
#build matchdf
nmat <- nrow(matches)
# make starter df
seen<-1*t((matches>0)[,2:3] )
mstat<-cbind(c(seen),c(1,2))
for(x in 1:nmat) {
if (matches[x,1] ==1) mat<-rbind(Pdf[matches[x,2],],Sdf[matches[x,3],])
 else if (matches[x,1] ==2) mat<-rbind(Pdf[matches[x,2],],Pdf[matches[x,2],])
else   mat<- rbind(Sdf[matches[x,3],],Sdf[matches[x,3],])
if (x==1)   matdf<-mat
else matdf<-rbind(matdf,mat)
}  #end for
matchdf <- cbind(mstat,matdf )
names(matchdf)[1] <- "seen"
names(matchdf)[2] <- "PorS"
return(matchdf)
}
