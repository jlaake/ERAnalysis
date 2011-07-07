\name{EarlyEffort}
\alias{EarlyEffort}
\docType{data}
\title{Effort Data for Gray Whale Southbound Migration Survey Data from 1967-1985}
\description{
Daily log of effort during the southbound migration survey of gray whales from 1967-1985.}
\usage{data(EarlyEffort)}
\format{
  A data frame with 1765 records on the following 23 variables.
  \describe{
    \item{\code{ID}}{Auto-number from Access database}
    \item{\code{Location}}{Location code: 1 is Yankee Point and 2 is Granite Canyon; different from Location in ERSurveyData}
    \item{\code{Day}}{Numeric day of month}
    \item{\code{Month}}{Numeric month}
    \item{\code{Year}}{4 digit year}
    \item{\code{Watchperiod}}{Watch period: 1,2 or 3; usually 1:0700-1200, 2:1200-1700; rarely 3}
    \item{\code{Observer}}{Observer number - see \code{\link{Observer}} table}
    \item{\code{Begintime}}{Beginning time of watch in HHMM format}
    \item{\code{Endtime}}{Ending time of watch in HHMM format}
    \item{\code{Block1winddirection}}{First recorded WindDirection (0 to 360)}
    \item{\code{Block1beaufort}}{First recorded Beaufort}
    \item{\code{Block1visibility}}{First recorded visibility code 1: Ex, 2: Very Good, 3: Good, 4: Fair, 5:Poor, 6: Unusable}
    \item{\code{Block2winddirection}}{Second recorded WindDirection (0 to 360)}
    \item{\code{Block2beaufort}}{Second recorded Beaufort}
    \item{\code{Block2visibility}}{Second recorded visibility code 1: Ex, 2: Very Good, 3: Good, 4: Fair, 5:Poor, 6: Unusable}
    \item{\code{Block3winddirection}}{Third recorded WindDirection (0 to 360)}
    \item{\code{Block3beaufort}}{Third recorded Beaufort}
    \item{\code{Block3visibility}}{Third recorded visibility code 1: Ex, 2: Very Good, 3: Good, 4: Fair, 5:Poor, 6: Unusable}
    \item{\code{key}}{key field in format yyyy-mm-dd_w where w is watch period; used to match effort with sightings}
    \item{\code{Use}}{set to "Yes" if effort/sightings all have a vis and beaufort that match specified maximums which has been 4.}
    \item{\code{Begin.date.time}}{Date-time of beginning of watch in POSIXct format}
    \item{\code{End.date.time}}{Date-time of end of watch in POSIXct format}
    \item{\code{Start.year}}{Numeric 4 digit year at the beginning of the survey}
  }
}
\details{
These data are from the EarlyEffort table in the Access database for gray whales.  They are effort records for the early surveys
from 1967 to 1985 which were conducted with a single effort record per watch so there is no clear way to split the effort into periods of 
common environmental values (eg vis, beaufort). These data are re-formatted into a common format and combined with the more recent survey data
into the \code{\link{PrimaryEffort}} dataframe.}
\keyword{datasets}
