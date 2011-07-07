\name{psfit.gamma}
\Rdversion{1.2}
\alias{psfit.gamma}
\alias{gammad}
\title{Computes probabilities and negative log-likelihood for truncated and discretized gamma distribution over range 0 to nmax
where range x-1 to x is mapped to integer x value for pod size.}
\description{
The \code{gammad} function provides probabilities for the truncated discretized gamma distribution
and the \code{psfit.gamma} function computes the negative log-likelihood for
the parameters of the gamma distribution for the observed
pod sizes for a given true pod size or range of true pod sizes from the calibration
data.}
\usage{
psfit.gamma(par, psmat, True = NULL, shape = TRUE, nmax=20)
gammad(par, nmax, True = NULL, shape = TRUE)
}
\arguments{
  \item{par}{log of gamma parameters; see details}
  \item{nmax}{maximum pod size}
  \item{psmat}{matrix of observed pod size frequencies in calibration experiment data;
               one row per pod of a known size and the value cell is the number of
               recorded pods of that size where the columns range from recorded size 1 to nmax}
  \item{True}{vector of true pod size for a range of true pod sizes when fitting a
              plus group (eg pods of size 4 and greater) which must match
              number of rows in \code{psmat}}
  \item{shape}{if TRUE then for plus groups the relationship with True pod size
              is for shape parameter and if FALSE it is used with scale=1/rate
              parameter}
}
\details{
If \code{True} is NULL then log(shape)=par[1], log(rate)=par[2] and scale=1/rate.
If \code{True} is not null and !shape then log(shape)=par[1], log(rate)=-par[2]-par[3]*True.
If \code{True} is not null and shape is TRUE then log(shape)=par[2]+par[3]*True, log(rate)=par[1].

Note that \code{psfit.gamma} is no longer used.  Instead refer to \code{\link{fit.gamma}}.
}
\value{
\code{gammad} returns a vector of probabilites for each size from 1 to nmax
\code{psfit.gamma} returns negative log-likelihood for values of parameters
for distribution of recorded pod sizes. }

\author{Jeff Laake}
\seealso{\code{\link{psfit.poisson}},\code{\link{psfit.negbin}}}

