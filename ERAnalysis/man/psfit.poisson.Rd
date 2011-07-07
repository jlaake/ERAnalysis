\name{psfit.poisson}
\name{psfit.poisson}
\Rdversion{1.1}
\alias{psfit.poisson}
\alias{poisson.d}
\title{Fitting code for truncated Poisson distribution with pod size calibration data}
\description{code{psfit.poisson} is no longer used but \code{poisson.d} is still useful to
compute probabilities for the truncated Poisson distribution.}
\usage{
psfit.poisson(par, psmat, True = NULL, nmax=20)
poisson.d(par, nmax, True = NULL)
}
\arguments{
  \item{par}{if True is Null then par[1]=log(lambda) where lambda is the Poisson intensity parameter;
  otherwise log(lambda)=par[1]+par[2]*True }
  \item{nmax}{maximum pod size}
  \item{psmat}{matrix of observed pod size frequencies in calibration experiment data;
               one row per pod of a known size and the value cell is the number of
               recorded pods of that size where the columns range from recorded size 1 to nmax}
  \item{True}{vector of true pod size for a range of true pod sizes when fitting a
              plus group (eg pods of size 4 and greater) which must match
              number of rows in \code{psmat}}
}
\details{
\code{poisson.d} was used by the deprecated function \code{psfit.poisson} and does nearly
the same thing as \code{\link{pois.d}} but the latter is more flexible.

Deprecated function: \code{psfit.poisson} function computes the negative log-likelihood for
the parameters of the truncated Poisson distribution for the observed
pod sizes for a given true pod size or range of true pod sizes from the calibration
data.
}
\value{
\code{poisson.d} returns a vector of probabilites for each size from 1 to nmax
\code{psfit.poisson} returns negative log-likelihood for values of parameters
for distribution of recorded pod sizes. a vector of probabilites for each size from 1 to nmax}

\author{Jeff Laake}
