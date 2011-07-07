\name{gsS}
\alias{gsS}
\alias{add.cf.all}
\alias{add.cf.laake}
\alias{add.cf.reilly}
\docType{data}
\title{Pod Size Correction Statistics}
\description{
Summarized calibration data for pod size estimation error. \code{gsS} is a calibration matrix
and the others are additive estimates.}
\usage{
data(gsS)
data(add.cf.all)
data(add.cf.reilly)
data(add.cf.laake)
}
\format{             
  \describe{
    \item{\code{gsS}}{A matrix with 20 rows and columns; each row is the true size and
                      each column the estimated size. Value is probability that a pod
                      of a true size will be estimated to be a particular size.}
    \item{\code{add.cf.all}}{vector of four additive correction factors (size 1,2,3 and 4+) using all of the
                             pod size calibration data via the Reilly approach.}
    \item{\code{add.cf.reilly}}{vector of four additive correction factors (size 1,2,3 and 4+) using 1978/79
                             pod size calibration data via the Reilly approach.}
    \item{\code{add.cf.all}}{vector of four additive correction factors (size 1,2,3 and 4+) using 1992/93 and 1993/94
                             pod size calibration data via the Reilly approach.}
    }
}
\details{
See \code{\link{create.podsize.calibration.matrix}} and \code{\link{reilly.cf}} for details on calculation for the values contained
within these computed data sets. The number of rows and columns in \code{gsS} depends on the
value set for \code{nmax}, the maximum possible true/estimated pod size.  The default value is 20.
}
\keyword{datasets}
