\name{PrimaryOff}
\alias{PrimaryOff}
\docType{data}
\title{Primary Observer Off-effort Southbound Sightings Data for Gray Whale Southbound Migration Survey Data from 1987-2006}
\description{
Data for sightings of off-effort souhtbound gray whale pods by primary observer during the southbound migration surveys of gray whales from 1987-2006.
These are sightings  in which the t241 occurred outside the watch period.}
\usage{data(PrimaryOff)}
\format{
  A data frame with 321 records on the following 17 variables.
  \describe{
    \item{\code{day}}{Numeric day of month}
    \item{\code{month}}{Numeric month}
    \item{\code{year}}{4 digit year}
    \item{\code{etime}}{Event decimal time}
    \item{\code{watch}}{Not on watch so value is -9}
    \item{\code{t241}}{Computed time (decimal days since 1 Dec) at which whale crosses line extending from shed at angle of 241; Used south or north time prior 1979; 1979+ caculated value using 3.24 nm whale speed}
    \item{\code{distance}}{Distance offshore in kilometers}
    \item{\code{podsize}}{Recorded number of whales in pod}
    \item{\code{observer}}{All NA}
    \item{\code{vis}}{All NA}
    \item{\code{beaufort}}{All NA}
    \item{\code{wind.direction}}{All NA}
    \item{\code{effort.period}}{All NA}
    \item{\code{pods.per.hour}}{All NA}
    \item{\code{station}}{All 'P'}
    \item{\code{seen}}{All 1}   
    \item{\code{Start.year}}{Numeric 4 digit year at the beginning of the survey}
  }
}
\keyword{datasets}
