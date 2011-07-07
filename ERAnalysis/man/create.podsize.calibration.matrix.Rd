\name{create.podsize.calibration.matrix}
\Rdversion{1.2}
\alias{create.podsize.calibration.matrix}
\alias{create.gsS}
\alias{reilly.cf}
\alias{expected.podsize}
\alias{podsize.computations}
\title{Various functions for pod size error/calibration}
\description{
\code{create.podsize.calibration.matrix} creates \code{\link{gsS}} the pod size calibration matrix from the 
calibration data (\code{\link{PodsizeCalibrationData}}) and optionally saves to package directory.  Previous
versions did the fitting as well but this is now done separately with the addition of random effects models.
This function is now just a wrapper that calls \code{create.gsS} with the results stored in \code{gamma.pod} to
compute the values of gsS. It also computes the additive correction factor 
vectors \code{\link{add.cf.all}}, \code{\link{add.cf.reilly}}, and \code{\link{add.cf.laake}} and
saves all 4 data sets in the package data directory.  The additive correction vectors can be used with \code{reilly.cf} to create pod size 'corrections'
using the approach first suggested by Reilly (1981). With \code{expected.podsize} the true pod size distribution \code{fS} and the expected pod size can be computed from the gamma distribution
specified by the parameter argument are computed as well as the conditional distribution of the true size \code{S} given an observed size \code{s} 
(\code{fSs}) and conditional expectation \code{ESs} based on the pod size calibration matrix \code{gsS}.  The function \code{podsize.computations} 
creates some plots of observed pod size distribution and fits a gamma distribution for the true pod size distribution for each year from the input
pod size calibration matrix \code{\link{gsS}} and the distribution of observed pod sizes.
}
\usage{
create.podsize.calibration.matrix(package.dir="C:/Users/Jeff Laake/Desktop/MyDocuments/R Development/ERAnalysis/ERAnalysis/data",save=TRUE)
create.gsS(ps.results,nmax=20,True=NULL)
reilly.cf(x,add.cf)
expected.podsize(spar,gsS,nmax=20)
podsize.computations(x,gsS,nmax=20,plot=FALSE)
}
\arguments{
  \item{package.dir}{Path to the data directory of the \code{ERAnalysis} package}
  \item{save}{if TRUE writes out dataframes to package data directory}  
  \item{ps.results}{result list from \code{\link{fit.gamma.re}} for fitted gamma with pod random effect}
  \item{True}{vector of true pod size values}
  \item{x}{Vector of observed pod sizes for \code{reilly.cf} and dataframe of observations for \code{podsize.computations}}
  \item{add.cf}{Vector of additive correction factors for observed size 1,2,3,4+; see \code{\link{add.cf.all}} as an example}
  \item{spar}{Vector of 2 parameters for a gamma distribution representing the true pod size distribution}
  \item{nmax}{Maximum value for true value/estimate of a pod size}
  \item{gsS}{Pod size calibration matrix; see \code{\link{gsS}}}
  \item{plot}{if TRUE produces plots of fitted and observed pod size distributions}
}
\details{
For details on fitting functions see eqns 6-11 in Laake et al.(2009) and refer to functions:
\code{\link{fit.gamma.re}},\code{\link{fit.gamma}},\code{\link{fit.poisson.re}},\code{\link{fit.poisson}}.
Previous functions \code{\link{psfit.gamma}},\code{\link{psfit.poisson}},\code{\link{psfit.negbin}} are
no longer used but have been left in the package.  Previously the calibration models were fit individually 
for a single categrory of true size (S=1,2,3,4+) with the data in tabular form (\code{\link{PodsizeCalibrationTable}}
but that code has been replaced with other functions that fit fixed or mixed (single random component) models
using all of the data (i.e., all values of S) simultaneously using \code{\link{PodsizeCalibrationData}}.

The additive correction factor approach is described in the Appendix of Laake et al (Additive Pod Size Correction Factor) and is
implemented with the function \code{reilly.cf} with the resulting vectors in vectors \code{\link{add.cf.all}}, \code{\link{add.cf.reilly}}, and \code{\link{add.cf.laake}}
for all of the data pooled, 1978/79 data and 1992/93 and 1993/94 respectively.
}
\value{
\code{create.podsize.calibration.matrix} does not return a value.

\code{reilly.cf} returns a vector of pod sizes that have been 'corrected' with the appropriate additive correction factor

\code{expected.podsize} returns a list with the following items:
    \item{\code{fS}}{Vector of probabilities for the true pod size distribution from 1 to \code{nmax}}
    \item{\code{ES}}{Value for the expected true pod size}
    \item{\code{fSs}}{Matrix of probabilities for the conditional true (row) pod size distribution from 1 to \code{nmax} given an observed value (column) from 1 to \code{nmax}}
    \item{\code{ESs}}{Vector of expected true pod sizes given observed pod sizes 1 to \code{nmax}}

\code{podsize.computations} returns a list of fitted gamma distributions for true pod size using \code{\link{fit.fS}}
which assumes pod size does not affect detectability.  These fitted distributions are used for starting values for
the analysis that does allow for pod size to affect detection (\code{\link{fit.pods}} and \code{\link{fit.missed.pods}}). 

}
\author{Jeff Laake}
\seealso{\code{\link{fit.gamma.re}},\code{\link{fit.gamma}},\code{\link{fit.poisson.re}},\code{\link{fit.poisson}}}
\references{LAAKE, J., A. PUNT, R. HOBBS, M. FERGUSON, D. RUGH, and J. BREIWICK. 2009. Re-evaluation of gray whale southbound surveys. NOAA Tech. Memo.
}
 