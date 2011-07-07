\name{psfit.negbin}
\Rdversion{1.1}
\alias{psfit.negbin}
\alias{negbin.d}
\title{Computes probabilities and negative log-likelihood for truncated negative binomial
range 1 to nmax.}
\description{
The \code{negbin.d} function provides probabilities for the truncated negative binomial distribution
and the \code{psfit.negbin} function computes the negative log-likelihood for
the parameters of the negative binomial distribution for the observed
pod sizes for a given true pod size or range of true pod sizes from the calibration
data.}
\usage{
negbin.d(par, nmax, True = NULL)
psfit.negbin(par, psmat, True = NULL, nmax=20)
}
\arguments{
  \item{par}{negative binomial parameters; see details}
  \item{nmax}{maximum pod size}
  \item{psmat}{matrix of observed pod size frequencies in calibration experiment data;
               one row per pod of a known size and the value cell is the number of
               recorded pods of that size where the columns range from recorded size 1 to nmax}
  \item{True}{vector of true pod size for a range of true pod sizes when fitting a
              plus group (eg pods of size 4 and greater) which must match
              number of rows in \code{psmat}}
}
\details{
If \code{True} is NULL then log(size)=par[1], logit(p)=par[2].
If \code{True} is not NULL then log(size)=par[1]+par[2]*True, logit(p)=par[3].
}

\value{
\code{negbin.d} returns a vector of probabilites for each size from 1 to nmax
\code{psfit.negbin} returns negative log-likelihood for values of parameters
for distribution of recorded pod sizes. }
\author{Jeff Laake}
\seealso{\code{\link{psfit.poisson}},\code{\link{psfit.gamma}}}

