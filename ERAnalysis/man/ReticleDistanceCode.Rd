\name{RetDistK}
\Rdversion{1.1}
\alias{RetDistK}
\alias{D241}
\alias{T241H}
\title{Compute Distance from Reticle and Offshore Distance and Time Abeam}
\description{
Functions to compute distance from reticles and offshore perpendicular distance from
bearing and distance and time to cross abeam of the shore based counting station. }
\usage{
RetDistK(Height, RadPerReticle=0.00497, Reticles)  
D241(sangle, sdist) 
T241H(stime, sangle, sdist , Speed=6.0)  
}
\arguments{
  \item{Height}{Height of observation platform in meters}
  \item{RadPerReticle}{Radians per reticle}
  \item{Reticles}{Number of reticles for measurement}
  \item{sangle}{Magnetic bearing in degrees of sighting of gray whales}
  \item{sdist}{Distance measurement from observer to whale at time of sighting}
  \item{stime}{Decimal time of sighting}
  \item{Speed}{Whale travel speed in km/hr}
}
\value{
\code{RetDistK} returns a vector of distances in kilometers.

\code{D241} returns vector of offshore perpendicular distances.

\code{T241H} returns vector of decimal times at which whales would be abeam (241 bearing)
of the shore based counting station.
}
\author{Rod Hobbs, Jeff Laake}
