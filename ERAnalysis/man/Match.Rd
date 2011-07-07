\name{Match}
\alias{Match}
\docType{data}
\title{Matched Sighting Data for Periods of Double-Observer Gray Whale Southbound Migration Survey Data from 1987-2006}
\description{
Matched sighting data for primary and secondary observer when more than one observer.  The data contain a record for the primary and secondary
observer for each uniquely defined pod that was seen by at least one of the observers. These data are used to construct a detection probability to 
predict the number of pods missed by one or both observers. The number of records can vary depending on the matching criterion used to construct the data.
Values of recorded sighting data can vary by \code{station} when bother detect a pod.}
\usage{data(Match)}
\format{
  A data frame with 22432 records on the following 19 variables.
  \describe{
    \item{\code{Date}}{yyyy-mm-dd character string}
    \item{\code{seen}}{1 if seen and 0 if missed}
    \item{\code{station}}{P if primary observer and S if secondary observer}
    \item{\code{day}}{Numeric day counting from 1 December of beginning of survey}
    \item{\code{watch}}{Watch period (1,2,3) based on time abeam of survey location (\code{t241})}
    \item{\code{t241}}{Computed time at which whale crosses line extending from shed at angle of 241; Used south or north time prior 1979; 1979+ caculated value using 3.25 nm whale speed}
    \item{\code{distance}}{Distance offshore in kilometers}
    \item{\code{podsize}}{Recorded number of whales in pod}
    \item{\code{Observer}}{observer code or initials; not complete for all years without matching to effort period; see \code{\link{Observer}} table}
    \item{\code{vis}}{At time of sightings visibility code 1: Ex, 2: Very Good, 3: Good, 4: Fair, 5:Poor, 6: Unusable}
    \item{\code{beaufort}}{Beaufort code at time of sightings}
    \item{\code{wind.direction}}{Wind direction at time of sighting as either a bearing (0-360) or compass directions (eg NW)}
    \item{\code{pphr}}{Number of observed pods passing per hour by the primary station during the effort period containing the sighting}
    \item{\code{mscore}}{Total matching score value; -99 if not matched}
    \item{\code{mcode}}{Not sure what this is - from ERAbund program}
    \item{\code{Start.year}}{Numeric 4 digit year at the beginning of the survey}
    \item{\code{original.watch}}{Original watch period in data; set presumably by event time for early surveys}
    \item{\code{hours}}{Number of hours of watch experience prior to current date}  
    \item{\code{Sex}}{Sex of the observer (M/F)}          
  }
}
\details{
Note that while the dataframe \code{\link{Match}} created by ERAbund is still maintained in the
package, those data are not used in the analysis because there was no way to automate changes in matching parameters and a few minor errors
were discovered in the ERAbund matching code (e.g., unless both observers started at the exact same time the data were not used). 
Instead, \code{\link{create.match}} and accompanying functions were created to build the match data.  It provides much more flexibility in
the type of criterion and in the weights used in the criterion.  
}
\keyword{datasets}
