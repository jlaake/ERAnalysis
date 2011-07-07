\name{SecondarySightings}
\alias{SecondarySightings}
\docType{data}
\title{Secondary Sighting Data for Gray Whale Southbound Migration Survey Data from 1987-2006}
\description{
Sighting data for secondary observer (when more than one observer) during the southbound 
migration surveys of gray whales from 1987-2006. 
}
\usage{data(SecondarySightings)}
\format{
  A data frame with 9726 records on the following 13 variables.
  \describe{
     \item{\code{day}}{Numeric day of month}
    \item{\code{month}}{Numeric month}
    \item{\code{year}}{4 digit year}
    \item{\code{etime}}{Event time in decimal days from 1 Dec of survey starting year}
    \item{\code{t241}}{Computed time at which whale crosses line extending from shed at angle of 241; Used south or north time prior 1979; 1979+ caculated value using 3.25 nm whale speed}
    \item{\code{distance}}{Distance offshore in kilometers}
    \item{\code{podsize}}{Recorded number of whales in pod}
    \item{\code{vis}}{At time of sightings visibility code 1: Ex, 2: Very Good, 3: Good, 4: Fair, 5:Poor, 6: Unusable}
    \item{\code{beaufort}}{Beaufort code at time of sightings}
    \item{\code{wind.direction}}{Wind direction at time of sighting as either a bearing (0-360) or compass directions (eg NW)}
    \item{\code{off}}{NA if on-watch; otherwise 'Off Watch'}
    \item{\code{Start.year}}{Numeric 4 digit year at the beginning of the survey}
    \item{\code{date}}{yyy-mm-dd date of sighting}
  }
}
\keyword{datasets}
