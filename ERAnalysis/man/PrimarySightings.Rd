\name{PrimarySightings}
\alias{PrimarySightings}
\docType{data}
\title{Primary Observer Filtered Southbound Sightings Data for Gray Whale Southbound Migration Survey Data from 1967-2006}
\description{
Data for sightings of souhtbound gray whale pods by primary observer during the southbound migration surveys of gray whales from 1967-2006. A common format is used to merge the data from the
early and recent surveys.  These are only sightings that occurred during the defined vis/beaufort conditions and in which the t241 occurred
during the watch period.}
\usage{data(PrimarySightings)}
\format{
  A data frame with 36244 records on the following 19 variables.
  \describe{
    \item{\code{Date}}{yyyy-mm-dd character string}
    \item{\code{day}}{Numeric day of month}
    \item{\code{month}}{Numeric month}
    \item{\code{year}}{4 digit year}
    \item{\code{watch}}{Watch period (1,2,3) based on time abeam of survey location (\code{t241}); this is based on more recent definitions which are 1: <10:30, 2:10:30-13:30, 3: >13:30 for all surveys even though earlier surveys only had 2 watches}
    \item{\code{t241}}{Computed time (decimal days since 1 Dec) at which whale crosses line extending from shed at angle of 241; Used south or north time prior 1979; 1979+ caculated value using 3.25 nm whale speed}
    \item{\code{distance}}{Distance offshore in kilometers}
    \item{\code{podsize}}{Recorded number of whales in pod}
    \item{\code{Observer}}{Observer code or initials; see \code{\link{Observer}} table}
    \item{\code{vis}}{Visibility code at time of sighting,  1: Ex, 2: Very Good, 3: Good, 4: Fair, 5:Poor, 6: Unusable}
    \item{\code{beaufort}}{Beaufort code at time of sighting}
    \item{\code{wind.direction}}{Wind direction at time of sighting as either a bearing (0-360) or compass directions (eg NW)}
    \item{\code{key}}{Character key field that links sightings to effort; has different formats for early & recent surveys}
    \item{\code{pphr}}{Number of observed pods passing per hour during the watch period containing the sighting}
    \item{\code{Start.year}}{Numeric 4 digit year at the beginning of the survey}
    \item{\code{original.watch}}{Original watch period in data; set presumably by event time for early surveys}
    \item{\code{only}}{TRUE if observation made when only Primary observer on watch} 
    \item{\code{hours}}{Number of hours of watch experience prior to current date}        
    \item{\code{Sex}}{Sex of the observer (M/F)}        
  }
}
\details{
Below is the table of sightings that match the filter of vis<=4 and beaufort<=4 which was
used by Buckland and Breiwick (2002) for the analysis of the 1967-1987 data.
Their totals for the years from 1967-1985 in Table 3 were: 903,1079,1245,1458,857,1539,1496,1508,1187,1991,657,
1730,1451,1756,1796. Some of the discrepancies are due to edits of the data that were found subsequent to that
publication. However the largest discrepancy of 44 sightings is for the 1985 survey when there were no edits made.
This is likely due to the restriction for beaufort and visibility for an entire watch.  For the data
used here the watch is assigned to a sighting based on the time the sighting crossed the line abeam of
the counting station (t241). No description of a calculation of t241 is given in the Buckland and Breiwick (1993) paper.
For 1984/85 there were 84 modifications to watches based on t241 and for 1985/86 there were 94 which could easily explain the
differences.

\preformatted{
1967: 903     1979: 1457
1968: 1072    1984: 1736
1969: 1236    1985: 1840
1970: 1463    1987: 2817
1971: 859     1992: 1544
1972: 1539    1993: 2348
1973: 1497    1995: 1909
1974: 1508    1997: 2196
1975: 1188    2000: 1544
1976: 1992    2001: 1598
1977: 657     2006: 1606
1978: 1735
}

The Buckland et al. (1993) reported 2981 North sightings and 2988 South sightings which matches
the counts of sightings in ERSurveyData for EFLAG=3 (sighting) for that year and value of LOCATION.  The smaller
count below for the south station (primary in 1987) of 2817 excludes northbound whales and any seen in
effort periods with vis>4 or beaufort>4.  The smaller numbers of 2426 and 2404 for N and S in Buckland and Breiwick (2002)
are presumably due to application of the all/or nothing exclusion based on an entire watch period although that is
not clearly documented in the paper.  If you run the following code:

\preformatted{
xx=merge(PrimarySightings[PrimarySightings$Start.year==1987,],PrimaryEffort[PrimaryEffort$Start.year==1987,],by="key")
dim(xx[xx$Use,]) 
}

It will report 2370 sightings which is 34 less than the total reported by Buckland and Breiwick (2002) for the south station, but as 
stated earlier our assignment for the watch is based on the t241 and not the initial time of sighting.

The discrepancies with other reports for the more recent surveys is primarily due to use of a different
visibility and beaufort filter.  For example, Laake et al. (1993) used visibility <=3 and no filter on
Beaufort.  And subsequent analyses have used a filter of visibility <=4 and no beaufort filter.

}
\keyword{datasets}
