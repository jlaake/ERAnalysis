\name{ERSurveyData}
\alias{ERSurveyData}
\docType{data}
\title{Gray Whale Southbound Migration Survey Data from 1987-2006}
\description{
Daily log of effort and sightings during the southbound migration survey of gray whales at
Granite Canyon California from 1987-2006.}
\usage{data(ERSurveyData)}
\format{
  A data frame with 44533 records on the following 28 variables.
  \describe{
    \item{\code{SEQUENCE}}{Year-sequence number in the order they were contained in the original raw data file; not always chronological}
    \item{\code{LOCATION}}{Location at Granite Canyon (GC): S- South, N-North, B- Bigeye shed, P-GC parking lot; or at Point Lobos: L,LS,LN}
    \item{\code{EXPERIMENT}}{Code: 1 - Single observer survey, 2 - Dual observer survey,  3 - Linkage test  4 - Distant whale, 5 - Pod size calibration, 6 - Whale tracking}
    \item{\code{DATE}}{Date}
    \item{\code{OBSERVER}}{Observer initials or code - see \code{\link{Observer}} table}
    \item{\code{EFLAG}}{Event flag: 1 - Start watch, 2-Environment change, 3 - gray whale sighting, 4 - comment, 5 - End watch, 6 - Other cetacean sighting, 7 - Start vessel watch, 8 - Vessel sighting, 9 - End vessel watch}
    \item{\code{ETIME}}{Time - decimal hours}
    \item{\code{VISCODE}}{Visibility code 1: Ex, 2: Very Good, 3: Good, 4: Fair, 5:Poor, 6: Unusable}
    \item{\code{WINDDIR}}{WindDirection (eg NW, SSW etc or 0-360 )}
    \item{\code{WINDFORCE}}{Beaufort code}
    \item{\code{NTIME}}{Time of north sighting (decimal hours)}
    \item{\code{NANGLE}}{Magnetic angle of north sighting}
    \item{\code{NRET}}{North sighting reticle measurement using 7x50 power binoculars; prior to 1995 only one reticle measurement stored in NRET and SRET}
    \item{\code{NEST}}{E if estimated reticle; For 1987/88 estimated distances were converted to reticle value}
    \item{\code{STIME}}{Time of south sighting (decimal hours)}
    \item{\code{SANGLE}}{Magnetic angle of south sighting}
    \item{\code{SRET}}{South sighting reticle measurement using 7x50 power binoculars; prior to 1995 only one reticle measurement stored in NRET and SRET}
    \item{\code{SEST}}{E if estimated reticle; For 1987/88 estimated distances were converted to reticle value}
    \item{\code{PODSIZE}}{Number of animals in pod}
    \item{\code{C_CPAIR}}{Number of cow/calf pairs}
    \item{\code{TRAVELDIR}}{Whale travel direction: S - south, N - north, M - milling}
    \item{\code{BEHAVIOR}}{Behavior code: 1 - breaching, 2 - spyhop, 3 - lobtailing, 4 - flipper  slap, 5 - mating, 6 - feeding, 7 - resting, 8 - associated with other animals, 9 - porpoising}
    \item{\code{NCUES}}{Number of cues}
    \item{\code{COMMENTS}}{Observer comments; in 1987/88 comments constructed with watch, distance class if estimated and distance value used if reticle}
    \item{\code{Start.year}}{Numeric 4 digit year at the beginning of the survey}
    \item{\code{year}}{Year (4 digit) of the observation}
    \item{\code{month}}{Numeric month of the observation}
    \item{\code{day}}{Numeric day of the month of observation}
  }
}
\details{
These data are from the AllRecentData table in the Access database for gray whales.  They are for the recent surveys
from 1987 to the present which all use the same format of a daily log with events to begin/end effort, sightings etc. These
data are used to create the \code{\link{PrimarySightings}} and \code{\link{PrimaryEffort}} dataframes which combine the effort and 
sightings for all of the gray whales surveys from 1967 to the present.
}
\keyword{datasets}
