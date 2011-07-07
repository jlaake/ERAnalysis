\name{fit.pods}
\Rdversion{1.1}
\alias{lnl.pods}
\alias{fit.pods}
\title{Estimates parameters for true pod size distribution}
\description{
Fits a gamma distribution for true pod size from records of a single observer
while accounting for missed pods and pod size error with a detection function
and pod size calibration matrix estimated from external data.}
\usage{
lnl.pods(par, formula, dpar, data, gsS, debug)
fit.pods(par, formula, dpar, data, gsS, debug, hessian)
}
\arguments{
  \item{par}{vector of 2 parameter values for pod size gamma distribution; log(shape) and log(rate)}
  \item{formula}{the formula for the logistic detection model}
  \item{dpar}{detection function parameter values}
  \item{data}{records for observations of primary station}
  \item{gsS}{nmax x nmax pod size calibration matrix; each row is a true pod size
             from 1 to nmax and the value for each column is the probability that
             a pod of a true size S is recorded as a size s (1..nmax columns)}
  \item{debug}{if TRUE will show iteration values}
  \item{hessian}{if TRUE will return hessian in model output from optim}
}
\details{
\code{fit.pods} calls \code{optim} to minimize the negative log-likelihood for
the parameters of the true pod size distribution which is computed by \code{lnl.pods}
based on the observed data and the externally derived detection function and the
pod size calibration matrix.
}
\value{
  \item{lnl.pods}{returns negative log-likelihood value}
  \item{fit.pods}{returns optim list of fitting results}
}
\author{Jeff Laake}
\seealso{\code{\link{fit.fS}}}
