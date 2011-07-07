\name{SecondaryEffort}
\alias{SecondaryEffort}
\alias{AllSecondaryEffort}
\docType{data}
\title{Effort Data of Secondary observer for Gray Whale Southbound Migration Survey Data from 1967-2006}
\description{
Effort data for secondary observer (when more than one observer) during the southbound migration surveys of gray whales from 1987-2006. 
These are only effort records with vis <4 and beaufort<=4 except that AllSecondaryEffort uses beaufort<=10 and vis<=10 because this latter
dataframe is used by \code{\link{create.match}} to include all sightings regardless of the condition.}
\usage{
data(SecondaryEffort)
data(AllSecondaryEffort)
}
\format{
  A data frame with 1600 records (2124 for AllSecondaryEffort) on the following 15 variables.
  \describe{
    \item{\code{watch.key}}{yyyy-ddw key field, where yyyy is 4 digit Start.year, dd is days since 1 Dec and w is the watch}
    \item{\code{Start.year}}{Numeric 4 digit year at the beginning of the survey}
    \item{\code{key}}{Effort period key which links effort period to sightings}
    \item{\code{begin}}{Begin time of effort period in decimal days from 1 Dec of survey \code{Start.year} 4 digit year}
    \item{\code{end}}{End time of effort period in decimal days from 1 Dec of survey \code{Start.year} 4 digit year}
    \item{\code{npods}}{Number of pods (estimate or observed) during the effort period}
    \item{\code{nwhales}}{Number of whales (estimate or observed) during the effort period}
    \item{\code{effort}}{Length of effort period in decimal days}
    \item{\code{vis}}{Visibility code during watch: 1: Ex, 2: Very Good, 3: Good, 4: Fair, 5:Poor, 6: Unusable; for <1987 it is the average value}
    \item{\code{beaufort}}{Beaufort code during watch; for <1987 it is the average value}
    \item{\code{Observer}}{Observer code or initials; see \code{\link{Observer}} table}
    \item{\code{time}}{Mid-point of effort period}
    \item{\code{watch}}{watch period number: 1 to 3}
    \item{\code{Use}}{True or False to indicate whether effort period is used in the analysis; it is NA where vis or beaufort is NA}
    \item{\code{Date}}{Date yyyy-mm-dd character field}
  }
}
\details{The gray whale survey data are organized into one to three watches depending on the year.  One observer is on effort during the entire watch;
when there are double observers each observer has a watch with a separate Location (N/S). Prior to 1987, effort was only recorded as an entire
watch period but starting in 1987 the watch period was broken into segments with constant environmental data like visibility and wind force (beaufort).
These smaller segments are called effort periods. So for surveys before 1987 there is a single effort period for each watch and for
1987 and after there can be several effort periods per watch.  Because the environmental data could vary within a watch period prior to 1987, the
value used here is the average of the vis/beaufort values from the sightings during that watch if there are any sightings.  If there are no 
sightings during the watch then the average of the 3 vis/beaufort measurements are used (see \code{\link{EarlyEffort}}). There is
no secondary effort prior to 1987 because there was only a single observer for 1967-1985 surveys.
}
\keyword{datasets}
