\name{fit.fS}
\Rdversion{1.1}
\alias{fit.fS}
\title{Log-Likelihood for fitting true pod size distribution}
\description{
Computes negative log-likelihood at given estimates parameters of gamma distribution for true pod size from counts of
recorded pod sizes and an estimated calibration matrix.}
\usage{
fit.fS(par, pstable, gsS, nmax = 20)
}
\arguments{
  \item{par}{log of gamma parameter estimates}
  \item{pstable}{counts of recorded pod sizes}
  \item{gsS}{ nmax x nmax pod size calibration matrix; each row is a true pod size
             from 1 to nmax and the value for each column is the probability that
             a pod of a true size S is recorded as a size s (1..nmax columns)}
  \item{nmax}{maximum possible pod size}
}
\details{
Computes negative log-likelihood based on assumed multinomial distribution and
input pod size calibration matrix gsS.}
\value{Negative log-likelihood value which can be minimized via the optim function}
\author{
Jeff Laake
}
\seealso{\code{\link{lnl.pods}},\code{\link{podsize.computations}}
}
