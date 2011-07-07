##############################################################################
#  Function from geofunc.xla translated for R and changed to return Kilometers
# instead of NM
#############################################################################
RetDistK<- function (Height, RadPerReticle=0.00497, Reticles)  {
# Height in meters
# RaderReticle = Radians per reticle mark
# Reticles = number of reticles below horizon
# RetDist = distance in kilometers
    if(length(Height)==1)
      Height=rep(Height,length(Reticles))
    x <- sqrt(2 * 6366 * Height / 1000 + (Height / 1000) ^ 2)
    distance=x
    Angle=vector("numeric",length(x))
    whichNA=which(is.na(Reticles))
    Reticles[is.na(Reticles)]=0
    Angle[Reticles>0] <- atan(x[Reticles>0] / 6366)
    distance[Reticles>0] <- ((6366 + Height[Reticles>0] / 1000) * sin(Angle[Reticles>0] + 
         Reticles[Reticles>0]*RadPerReticle) - sqrt(6366 ^ 2 - ((6366 + Height[Reticles>0] / 1000) *
         cos(Angle[Reticles>0] + Reticles[Reticles>0] * RadPerReticle)) ^ 2))
    distance[whichNA]=NA
    return(distance)
}
##############################################################################
#  Function to calculate projected time to cross Azmuth at 241 degrees magnetic
#  from gray whale sighting uses swim speed, time at sighting, and recorded
# bino compass reticle and angle returns decimal hours
# note should add variance
# note speed and variance of speed are from  Swartz et al. 1986
#############################################################################
T241H <- function (stime, sangle, sdist , Speed=6.0)  {
    #stime  time of sighting in decimal hours
    #sangle magnetic angle from bino compass to sighting
    #dist  distance to sighting
    #Speed=6.0 km/hr
  return(stime + ((sin(pi*(sangle-241)/180)*sdist/Speed)))
}
###############################################################################
#  Function to calculate projected distance offshore at 241 degrees magnetic
#  from gray whale sighting uses recorded bino compass angle
# note should add variance
#############################################################################
D241 <- function (sangle, sdist)   {
    #sangle magnetic angle from bino compass to sighting
    #dist  distance to sighting
    #Speed=6.0 km/hr
 return(cos(pi*abs(sangle-241)/180)*sdist)
}