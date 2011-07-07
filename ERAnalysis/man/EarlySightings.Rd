\name{EarlySightings}
\alias{EarlySightings}
\docType{data}
\title{Sighting Data for Gray Whale Southbound Migration Survey Data from 1967-1985}
\description{
Daily log of sightings of gray whales during the southbound migration survey of gray whales from 1967-1985.}
\usage{data(EarlySightings)}
\format{
  A data frame with 28504 records on the following 32 variables.
  \describe{
    \item{\code{ID}}{Auto-number from Access database}
    \item{\code{Location}}{Location code: 1 is Yankee Point and 2 is Granite Canyon; different from Location in ERSurveyData}
    \item{\code{Day}}{Numeric day of month}
    \item{\code{Month}}{Numeric month}
    \item{\code{Year}}{4 digit year}
    \item{\code{Watchperiod}}{Watch period (1,2,3); usually 1:0700-1200, 2:1200-1700; rarely 3; originally assigned to sighting based on initial time of sighting}
    \item{\code{Northtime}}{Time of north sighting of gray whales in HHMM format}
    \item{\code{Southtime}}{Time of south sighting of gray whales in HHMM format}
    \item{\code{Podsize}}{Recorded number of whales in pod}
    \item{\code{Distanceoffshore}}{Estimated distance offshore in distance intervals; unclear whether it is statute or nautical miles- assumed to be nautical}
    \item{\code{Traveldirection}}{Direction of travel of the whales: S for south and N for North}
    \item{\code{Winddirection}}{Direction of the wind at time of sighting (0 to 360)}
    \item{\code{Beaufort}}{Beaufort code at time of sightings}
    \item{\code{Visibility}}{Visibility code at time of sighting: 1: Ex, 2: Very Good, 3: Good, 4: Fair, 5:Poor, 6: Unusable}
    \item{\code{Podsizeflag}}{* if no pod size in data and assumed value of 1 was entered}
    \item{\code{Northangle}}{Angle of north sighting (not recorded 1967-1977); in 1978 recorded +/- relative to abeam; 1979 on as magnetic bearing}
    \item{\code{Southangle}}{Angle of south sighting (not recorded 1967-1977); in 1978 recorded +/- relative to abeam; 1979 on as magnetic bearing}
    \item{\code{Verticalangle}}{Vertical angle ala reticle measurement at S. sighting}
    \item{\code{Behavior}}{Various behavioral codes used only in 1979}
    \item{\code{Timeinsert}}{Either N or S if north or south time was computed and inserted into the data record}
    \item{\code{Visflag}}{* if not in data and surmised from nearby sightings}
    \item{\code{Flukes}}{Number of flukes seen}
    \item{\code{Breaches}}{Number of breaches seen}
    \item{\code{Cue}}{Cue for sighting}
    \item{\code{Distanceflag}}{E - estimated; C- calculated}
    \item{\code{Sightabilitycode}}{1 or 2; reason for it not clear and only used in 1985/86 survey}
    \item{\code{Observer}}{Observer number; see \code{\link{Observer}} table}
    \item{\code{key}}{Key field of date and assigned watch period to link to effort}
    \item{\code{Assignedwatchperiod}}{Watch period assigned to period based on time abeam relative to begin/end times of watches}
    \item{\code{Included}}{Defined field which is True if Timeabeam is contained in a effort period and False otherwise}
    \item{\code{Timeabeam}}{Computed time (HH:MM:SS) at which whale crosses line extending from shed at angle of 241 abeam of 
                       shore station; Used south (if available) or north time prior to 1979 for Timeabeam; for 1979 and beyond
                       calculated value using 3.24 nm (6km/hr) whale speed}   
    \item{\code{Start.year}}{Numeric 4 digit year at the beginning of the survey}
  }
}
\details{
These data are from the EarlySightings table in the Access database for gray whales.  They are sighting records for the early surveys
from 1967 to 1985 which were conducted with a single effort record per watch so there is no clear way to split the effort into periods of 
common environmental values (eg vis, beaufort). These data are re-formatted into a common format and combined with the more recent survey data
into the \code{\link{PrimarySightings}} dataframe.}
\keyword{datasets}
