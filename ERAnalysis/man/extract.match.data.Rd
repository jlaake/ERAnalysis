\name{extract.match.data}
\Rdversion{1.1}
\alias{extract.match.data}
\title{Extract data for matching algorithm}
\description{
 Extracts the relevant data from recent surveys to do matching of double observer data.
 There are no arguments but it uses the dataframe ERSurveyData and Height which
 contains the platform heights at the observation locations each year. Excludes any
 observations of northbound whales.
}
\usage{
extract.match.data()
}
\value{
list of match data with each element being the Primary and Secondary observation of
gray whale sighings data for a particular day.
}
\author{Jeff Laake}


