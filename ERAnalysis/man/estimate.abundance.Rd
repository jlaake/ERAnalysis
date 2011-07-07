\name{estimate.abundance}
\Rdversion{1.1}
\alias{estimate.abundance}
\title{Horvitz-Thompson-like estimates of abundance for a migration survey}
\description{
Uses \code{\link{fit.migration.gam}} to fit a migration curve to a sequence of effort periods
and whale passage rates to create naive abundance estimates or abundance estimates
corrected for missed pods and pod size error.  The latter uses a fitted detection function 
and parameters for pod size distribution to construct estimates of abundance from the observations 
using a Horvitz-Thompson estimation approach (eg 1/p where p varies by observation). From the fitted
migration curve it computes total whale abundance.
}
\usage{
  estimate.abundance(spar,dpar,gsS,effort,sightings,dformula=~True,
                     gformula=~s(time),nmax=20,pod=FALSE,plotit=TRUE,anchor=TRUE,
                     show.anchor=FALSE,sp=NULL,final.time=90,lower.time=0,do.mult=TRUE,pool=TRUE,...)
}
\arguments{
  \item{spar}{ vector of 2 parameters for gamma distribution fitted to
               true pod size distribution}
  \item{dpar}{vector of parameters for detection function of pods as described
              by formula}
  \item{gsS}{nmax x nmax pod size calibration matrix; each row is a true pod size
             from 1 to nmax and the value for each column is the probability that
             a pod of a true size S is recorded as a size s (1..nmax columns)}
  \item{effort}{dataframe for observation effort; see details below}
  \item{sightings}{dataframe of observations for the primary observer which contains
              all of the fields (no NA's) used in formula; see details.}
  \item{dformula}{formula for the logit of detection probability created with the
                 double-observer data}
  \item{gformula}{formula for the quasi-Poisson GAM used to describe the migration curve over time}
  \item{nmax}{maximum possible pod size}
  \item{pod}{if TRUE fits to npods; otherwise, it uses nwhales in the dataframe}
  \item{plotit}{if TRUE will plot the fitted curve and observed data values}
  \item{anchor}{if TRUE will anchor curve at \code{lower.time} and \code{final.time} with values of 0}
  \item{show.anchor}{if TRUE and \code{anchor==TRUE & plot==TRUE} then non-anchored curve
                     is fitted and shown as a dashed line on the plot}
  \item{sp}{vector of smoothing parameters to be passed to gam; probably no longer needed with new version of mgcv}
  \item{final.time}{vector of final times to use for anchoring}  
  \item{lower.time}{vector of beginning times to use for anchoring}
  \item{do.mult}{if TRUE, the abundance is also computed based on a multiplier of the ratio of the
              total area under the migration curve divided by the area under the curve during the sampled
              periods}
  \item{pool}{if TRUE, all of the abundance data is pooled for the GAM fit which can be year-specific for some
              or all of the parameters}
  \item{...}{any additional arguments to be passed to gam}
}
\details{
The \code{sightings} and \code{effort} dataframes have some fairly specific requirements.
First and foremost, they must each have a field named \code{key} which is used to
link the \code{sightings} contained in a specific period of \code{effort}. As well,
\code{effort} must have fields named \code{Start.year} which is the year the
survey began, \code{begin} and \code{end} which are the
decimal day values for the begin and end times for the effort period,\code{time}
which is the mid-point of the interval and \code{effort} which is the length of
the interval. It can also have optional fields \code{vis} and \code{beaufort} which
are the values during the effort period.

The \code{sightings} dataframe must have a field \code{podsize}, the recorded observed
pod size, and optionally \code{corrected.podsize}. All other fields are optional but
there must be a valid value (not NA) for each field used in \code{dformula} for
computing detection probability.

The form of estimation depends on the values of \code{dpar}, \code{dformula} and
\code{spar}.  If all of these values are NULL, then there is no expansion for
missed pods or pod size correction using the calibration matrix \code{gsS}.  If
\code{pod=TRUE}, it fits the GAM to the observed number of pods in each effort period;
otherwise, it uses the sum of the \code{corrected.podsize} values (if available) per effort period
or it uses sum of the \code{podsize} values.  If \code{dpar} and \code{dformula} have
values then the number of pods or whales are corrected for detection probability using
a Horvitz-Thompson approach of 1/p_i or s/p_i where s is pod size and p_i is specific for
each of the observations indexed by i=1,...,n. If values are specified for \code{spar} they
are used for a gamma distribution for the true pod sizes. for each observation, the
conditional distribution of the true size given the observed size is computed and used
in the estimation.

Upon computing abundance for each effort period, these are combined with effort data
and used in a call to \code{\link{fit.migration.gam}} to estimate abundance throughout the
migration.

This function calls the following functions in the package: \code{\link{gammad}} and \code{\link{fit.migration.gam}}.
It contains the code for equations 17-19 in Laake et al.(2009).
}
\value{
  \item{result}{result list from \code{\link{fit.migration.gam}}}
  \item{options}{list of options used in call to this function}
}
\author{ Jeff Laake}


