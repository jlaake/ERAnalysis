\name{PodsizeCalibrationData}
\alias{PodsizeCalibrationData}
\alias{PSCalibAll}
\docType{data}
\title{Pod Size Calibration Data to Correct Pod Size Estimation Error}
\description{
Calibration data for pod size estimation error were recorded in 1978/79, 1993/93, and 1993/94 with aerial surveys and in
1997/98 with extended land-baed tracking.  One or more observers estimated the size of the pod which had a determined true size
from either the aerial observation or tracking.}
\usage{
data(PodsizeCalibrationData)
data(PSCalibAll)
}
\format{
  A data frame with 723 records on the following 7 variables.
  \describe{
    \item{\code{Type}}{Factor variable with value of either: Aerial or Tracking}
    \item{\code{Year}}{Factor variable with values: 1978/79, 1992/93, 1993/94, or 1997/1998}
    \item{\code{ID}}{ID assigned to pod within experiment}
    \item{\code{True}}{True size of pod}
    \item{\code{Observer}}{Observer initials for surveys other than 1978/79 which have arbitrary values that don't link to \code{\link{Observer}} table}
    \item{\code{Estimate}}{Observer's estimate of the pod size}
    \item{\code{key}}{Key field constructed of Type_Year_ID}
    }
}
\details{ These data are also in the text file PSCalibAll.txt in the data directory for the
package which is what is in PSCalibAll except for the key field.}
\keyword{datasets}
