\name{Height}
\alias{Height}
\docType{data}
\title{Eye heights of observers for surveys from 1967-2006}
\description{
Observer eye height (in meters) of observer for various locations for 1967-2006 for computation
of distance from reticle binoculars}
\usage{data(Height)}
\format{
  A data frame with 32 records on the following 4 variables.
  \describe{
    \item{\code{Start.year}}{4 digit year in which survey was started}
    \item{\code{Location}}{Location code: S- south, N- north, NP - north parking lot}
    \item{\code{key}}{Year or concatenated year and location value for merging observation data with heights; see \code{\link{extract.match.data}}}
    \item{\code{Height}}{Observer eye height in meters}
  }
}
\keyword{datasets}
