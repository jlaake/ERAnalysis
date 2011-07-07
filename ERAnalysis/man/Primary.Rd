\name{Primary}
\alias{Primary}
\docType{data}
\title{Primary Observer Southbound Sightings Data for Gray Whale Southbound Migration Survey Data from 1987-2006}
\description{
Data for sightings of souhtbound gray whale pods by primary observer during the 
southbound migration surveys of gray whales from 1987-2006.}
\usage{data(Primary)}
\format{
  A data frame with 17987 records on the following 16 variables.
  \describe{
    \item{\code{day}}{Numeric day of month}
    \item{\code{month}}{Numeric month}
    \item{\code{year}}{4 digit year}
    \item{\code{watch}}{Watch period (1,2,3) based on time abeam of survey location (\code{t241})}
    \item{\code{t241}}{Computed time (decimal days since 1 Dec) at which whale crosses line extending from shed at angle of 241; Used south or north time prior 1979; 1979+ caculated value using 3.24 nm whale speed}
    \item{\code{distance}}{Distance offshore in kilometers}
    \item{\code{podsize}}{Recorded number of whales in pod}
    \item{\code{Observer}}{Observer code or initials; see \code{\link{Observer}} table}
    \item{\code{vis}}{Visibility code at time of sighting,  1: Ex, 2: Very Good, 3: Good, 4: Fair, 5:Poor, 6: Unusable}
    \item{\code{beaufort}}{Beaufort code at time of sighting}
    \item{\code{wind.direction}}{Wind direction at time of sighting as either a bearing (0-360) or compass directions (eg NW)}
    \item{\code{key}}{Character key field that links sightings to effort; has different formats for early & recent surveys}
    \item{\code{pphr}}{Number of observed pods passing per hour during the watch period containing the sighting}
    \item{\code{Start.year}}{Numeric 4 digit year at the beginning of the survey}
    \item{\code{original.watch}}{Watch period (1,2,3) assigned by ERAbund based on t241}
    \item{\code{only}}{TRUE if observation made when only Primary observer on watch}    
  }
}
\keyword{datasets}
