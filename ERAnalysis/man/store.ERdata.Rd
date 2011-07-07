\name{store.ERdata}
\Rdversion{1.1}
\alias{store.ERdata}
\alias{store.ERAbundData}
\title{Retrieves and Stores Gray Whale Data into Package}
\description{
Retrieves data from Access database (\code{GrayWhaleSurveyData.accdb}) and from directories containing the \code{ERAbund} output files for the more recent surveys.  Using these 
dataframes it creates common files \code{\link{PrimarySightings}} and \code{\link{PrimaryEffort}} and others. These are saved into the 
data subdirectory for the package.  These functions are only needed if the data are changed or the criterion for data selection such as
the matching criterion or the visibility/beaufort criterion are changed.  After running these functions the package must be rebuilt to include
the new versions into the pacakge.}
\usage{
store.ERdata(DBDirectory='',package.dir='C:/Users/Jeff Laake/Desktop/MyDocuments/R Development/ERAnalysis/ERAnalysis/data',nmax=20)
store.ERAbundData(fdir='c:/gw',package.dir='C:/Users/Jeff Laake/Desktop/MyDocuments/R Development/ERAnalysis/ERAnalysis/data',
                                 maxbeauf=4,maxvis=4)
}
\arguments{
  \item{DBDirectory}{Path for the directory of the Access database; if blank it is assumed to be in the same directory as the current R workspace.}
  \item{package.dir}{Path to the data directory of the \code{ERAnalysis} package}
  \item{fdir}{Path to the main directory containing the sub-directories for the \code{ERAbund} data processing (eg a subdirectory would be named something like ERSW8788 where 8788 refers to the
               1987/88 survey}
  \item{maxbeauf}{Maximum value for Beaufort code that is allowed for effort and sightings merged into \code{\link{PrimarySightings}} and \code{\link{PrimaryEffort}}}
  \item{maxvis}{Maximum value for Visibility code that is allowed for effort and sightings merged into \code{\link{PrimarySightings}} and \code{\link{PrimaryEffort}}}
  \item{nmax}{Maximum value for true value/estimate of a pod size}
}  
\details{
\code{store.ERdata} retrieves data from the Gray Whale Access Database and creates dataframes
\code{\link{EarlyEffort}}, \code{\link{EarlySightings}} for surveys from 1967-1985 and \code{\link{ERSurveyData}} for
surveys from 1987-2006.  Some computed fields are added to the dataframes and then they are saved to the \code{ERAnalysis} package data directory.
The function also reads in the text file \code{PSCalibAll.txt} contained in the package data directory and creates the dataframes
\code{\link{PodsizeCalibrationData}} and  \code{\link{PodsizeCalibrationTable}} which are also saved in the package data directory.

\code{store.ERAbundData} calls \code{\link{process.ERdata}} to import and process the data output files from the ERAbund analysis processing
and merge the effort and sightings from the earlier surveys with the more recent surveys and saves them into the data directory for the package.
See \code{\link{process.ERdata}} for more details about the processing. If you make changes to the data in the Access database, then you 
need to use \code{store.ERdata} and rebuild the package before calling this function because it uses the \code{data} function to extract the
data from the package which would have been from the last package build.
}
\author{Jeff Laake}
\seealso{\code{\link{process.ERdata}},\code{\link{ReadERAbundFiles}}}
