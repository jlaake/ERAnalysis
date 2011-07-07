\name{ReadERAbundFiles}
\Rdversion{1.1}
\alias{ReadERAbundFiles}
\alias{process.ERdata}
\title{Reads files created by ERAbund data entry/analysis program}
\description{
\code{ReadERAbundFiles} creates dataframes Match, Primary, Secondary and ERNorm by reading in all
processed files from recent gray whale surveys and combining results across years. \code{process.ERdata} creates and converts some data and
merges the sightings and effort for all of the gray whale surveys from 1967-2006 that meet specified visibility and Beaufort criterion.
}
\usage{
ReadERAbundFiles(directory="c:/gw",years=c(87,92,93,95,97,0,1,6))
process.ERdata(fdir="c:/gw",maxbeauf=4,maxvis=4)
}
\arguments{
  \item{directory}{directory containing sub-directories of ERAbund result files}
  \item{years}{vector of last 2 digits of survey year; 87 is 1987/88 survey and 0 is 2000/2001 survey}
  \item{fdir}{Path to the main directory containing the sub-directories for the \code{ERAbund} data processing (eg a subdirectory would be named something like ERSW8788 where 8788 refers to the
               1987/88 survey}
  \item{maxbeauf}{Maximum value for Beaufort code that is allowed for effort and sightings merged into \code{\link{PrimarySightings}} and \code{\link{PrimaryEffort}}}
  \item{maxvis}{Maximum value for Visibility code that is allowed for effort and sightings merged into \code{\link{PrimarySightings}} and \code{\link{PrimaryEffort}}}
}
\details{
\code{process.ERdata} retrieves the survey data for 1967-1985 from the package data directory and
reads in the analysis output files from ERAbund using \code{ReadERAbundFiles}.  It creates various data fields in diffrent
formats and converts distance from meters to kilometers for the \code{ERAbund} files and from nautical miles to
kilometers for the early survey data (1967-1985).  It creates a common dataframe for migration curve data which
is eventually named \code{\link{PrimaryEffort}} by combining the effort data from \code{\link{EarlyEffort}} from the earlier surveys
and the effort data from the ERNorm files created by the ERAbund analysis from the more recent surveys. For the earlier surveys, the effort data is filtered such
that the effort for the entre watch is used as long as there was no recorded visibility or beaufort on the effort record or the sightings on that watch
that exceeded \code{maxbeauf} or \code{maxvis}.  If any value exceeds the set values, the entire watch is excluded.  The data field \code{Use} is 
defined based on the default values of 4 but if they are changed in the call to \code{process.ERdata} that value is over-written.
With the more recent data, effort periods are created within the watch that have a constant value for vis or beaufort, so instead of dropping entire
watch periods, the vis and beaufort criterion are applied to the much smaller effort period with their associated sightings. For the more recent
surveys, a field \code{Use} is defined and assigned a Yes/No value using the application of throwing out entire watch periods based on the criterion.
That will allow consistent application of the rule across all years if the analyst chooses.  But the effort periods for the more recent surveys
are included into \code{\link{PrimaryEffort}} based on the more restricted application of the exclusion criterion.  For the earlier surveys, there
is no single vis and beaufort values for the watch period, so averages are used of the values for the sightings during the watch period and if
there are no sightings then any effort block values are used (e.g., Block1visibility, Block2visibility, Block3visibility).  As with the effort,
the sightings from the primary observer are merged into a common format (\code{\link{PrimarySightings}}) for the early and more recent surveys and these are restricted to
effort periods that meet the vis and beaufort criterion. The field \code{pphr} which is pods per hour is computed by the ERAbund program for the
more recent surveys and is computed here for the early surveys using the observed number of pods/hour during the watch period in which the sighting
was seen.  In some cases, a distance was not recorded for a sighting during the earlier surveys.  As a default, the mean distance is used
of all sightings from the earlier surveys.
}
\value{
From \code{ReadERAbundFiles}
\item{Match}{dataframe of matched double-observer data from all processed years}
\item{Primary}{dataframe for all processed years of observation data from Primary observation station that met vis and beaufort
                criterion specified in the .in file for run of ERAbund analysis and was on-effort.}
\item{Secondary}{dataframe for all processed years of observation data from Secondary observation station that met vis and beaufort
                criterion specified in the .in file for run of ERAbund analysis.}
\item{ERNorm}{dataframe for all processed years of migration rate data that is set up
               in gwnorm structure.}

From \code{process.ERdata}
\item{Match}{dataframe of matched double-observer data from all processed years}
\item{Secondary}{dataframe for all processed years of observation data from Secondary observation station that met vis and beaufort
                criterion specified in the .in file for run of ERAbund analysis.}
\item{PrimarySightings}{dataframe of all sightings for primary observer for all survey years (1967-2006) that meet the
                        visibility and Beaufort criterion}
\item{PrimaryEffort}{dataframe of all effort for primary observer for all survey years (1967-2006) that meet the
                        visibility and Beaufort criterion}

}
\author{Jeff Laake}
\seealso{\code{\link{CreateERAbundFiles}},\code{\link{store.ERAbundData}}}
