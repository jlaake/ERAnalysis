\name{compute.sampling.multiplier}
\Rdversion{1.1}
\alias{compute.sampling.multiplier}
\alias{integrate.gam}
\title{Computes multiplier to expand from sampled effort to entire migration}
\description{
Integrates the area under the fitted GAM migration curve and the area under
the curve contained in the intervals that were sampled.  The ratio of the
total area to the area within the intervals is the multiplier to expand
the abundance for the entire migration. This implements equations 21-22 in Laake et al. (2009) }
\usage{
compute.sampling.multiplier(model,data,lower=0,upper=90)
integrate.gam(x,model,vis,beauf)
}
\arguments{
  \item{model}{fitted GAM model}
  \item{data}{migration curve dataframe}
  \item{lower}{beginning of migration curve;lower value for integration which is days from 1 Dec}
  \item{upper}{end of migration curve;upper value for integration which is days from 1 Dec}
  \item{x}{vector of times to be evaluated within range of migration curve}
  \item{vis}{vis values to use for prediction if needed}
  \item{beauf}{beauf values to use for prediction if needed}
}
\value{
\code{integrate.gam} returns a vector of response predictions from the GAM model that are
used with a call to integrate from \code{compute.sampling.multiplier}

\code{compute.sampling.multiplier} returns the multiplier which is the total
area under the migration curve divided by the area under the curve within the time intervals
that were sampled.
}
\author{Jeff Laake}

