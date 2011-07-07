\name{fit.migration.gam}
\Rdversion{1.1}
\alias{fit.migration.gam}
\title{Fit migration curve with a GAM}
\description{
Fits a GAM to the migration data over time using either number of pods or number
of whales.  It can optionally anchor the data to include 0 observations at time
\code{lower.time} and \code{final.time}.  It is called from \code{\link{estimate.abundance}} but can be 
called on its own.  It is used for estimation of naive abundance and for abundance corrected for missed pods and
pod size error.
}
\usage{
fit.migration.gam(er.migdata,years,formula=~s(time),pod=FALSE,plotit=TRUE,anchor=TRUE,
                                show.anchor=FALSE,sp=NULL,final.time=rep(90,length(years)),
                                lower.time=rep(0,length(years)),do.mult=TRUE, pool=TRUE,...)
}
\arguments{
  \item{er.migdata}{dataframe of migration data; see details}
  \item{years}{vector of survey years to be fitted}
  \item{formula}{formula for the GAM model}
  \item{pod}{if TRUE fits to npods; otherwise, it uses nwhales in the dataframe}
  \item{plotit}{if TRUE will plot the fitted curve and observed data values}
  \item{anchor}{if TRUE will anchor curve at \code{lower.time} and \code{final.time} with values of 0}
  \item{show.anchor}{if TRUE and \code{anchor==TRUE & plot==TRUE} then non-anchored curve
                     is fitted and shown as a dashed line on the plot}
  \item{sp}{vector of smoothing parameters to be passed to gam; probably no longer needed with new version of mgcv}
  \item{final.time}{vector of final times to use for anchoring}
  \item{lower.time}{vector of beginning times to use for anchoring}
  \item{do.mult}{if TRUE, computes multiplier to scale observed effort to full migration timeframe}
  \item{pool}{if TRUE, all of the abundance data is pooled for the GAM fit which can be year-specific for some
              or all of the parameters; if(pool) then do.mult is ignored.}
  \item{...}{any additional arguments to be passed to gam}
}
\details{
The dataframe must contain the following fields: 1) Start.year which is the 4 digit
year in which the survey started, 2) begin - the number of decimal days since 1 Dec for the beginning
of the effort period, 3) end - the number of decimal days since 1 Dec for the end of the effort period,
4) effort - length of period in decimal days, 5) either: npods - number of pods seen
or predicted to have passed during the effort period or nwhales - number of whales seen
or predicted to have passed during the effort period, 6) time - midpoint of the
period in decimal days since 1 Dec.

This function calls the following functions in the package: \code{\link{compute.sampling.multiplier}}. 
It contains the code for equations 25 and 26 in Laake et al.(2009).
}
\value{
\item{models}{list of fitted GAM models for each of the specified survey \code{years}}
\item{Total}{vector of estimated total abundances of either pods or whales computed as sum of
predictions for each day for each of the specified survey \code{years}}
\item{TotalM}{vector of estimated total abundance of either pods or whales computed as sum of
observations times the multiplier (described below) for each of the specified survey \code{years}}
\item{mult}{vector of multipliers for each of the specified \code{years}. The multiplier
is the ratio of the total area under the curve and the area under the curve within the sampled intervals
}
}
\author{
Jeff Laake
}
\seealso{\code{\link{compute.sampling.multiplier}}
}
