\name{fit.missed.pods}
\Rdversion{1.1}
\alias{fit.missed.pods}
\alias{lnl.missed.pods}
\title{Gray whale double-observer data fitting function}
\description{Fits true pod size distribution and the detection probability parameters
 for a given detection model specified by formula.}
\usage{
fit.missed.pods(data, primary, pbyyear=FALSE, formula = ~1,  par = NULL, 
          gsS, maxit = 1000, refit = TRUE, debug = FALSE, hessian = FALSE, method="BFGS")
lnl.missed.pods(par, data, primary.only, years=NULL, formula, gsS, debug)
}
\arguments{
  \item{data}{match data ordered in pairs with station=P then S where
              P=Primary and S=Secondary}
  \item{primary}{observations made by primary observer}
  \item{primary.only}{observations made by primary observer when there was no secondary observer}
  \item{pbyyear}{if TRUE, fits year-specific gamma parameters for pod size}
  \item{years}{vector of years included in the analysis}
  \item{formula}{the formula for the logistic detection model}
  \item{par}{initial parameter values including those for pod size distribution
              followed by those for detection model; if left blank initial estimates
              are derived - see details}
  \item{gsS}{nmax x nmax pod size calibration matrix; each row is a true pod size
             from 1 to nmax and the value for each column is the probability that
             a pod of a true size S is recorded as a size s (1..nmax columns)}
  \item{maxit}{maximum number of iterations}
  \item{refit}{if TRUE will continue to call optim to refit until convergence is achieved
               regardless of value of maxit}
  \item{debug}{if TRUE will show iteration values}
  \item{hessian}{if TRUE will return hessian in model output from optim}
  \item{method}{method to be used in optim for optimization}  
}
\details{The parameter vector (par) is ordered with the pod size distribution parameters
first followed by the parameters for detection probability. If \code{pbyyear=FALSE} there are
only 2 pod size parameters with the first for the shape and the second for the rate which
is 1/scale. The value of par is actually log(x) where x is the actual parameter in the
gamma distribution.  Thus, x=exp(par) which restricts the parameters of the gamma
to be strictly positive which is required. If \code{pbyyear=TRUE}, then there are 2*n
pod size parameters ordered in shape,rate pairs for each of the n ascending sorted years.
The detection function parameters follow the pod size parameters in the vector.  The
length of those parameters must match the number of colums in the design matrix
created by the \code{formula} with the \code{data}. In specifying initial parameters,
you can specify each value or only the values for pod size.  If initial values
are not provided then code is used to create starting values.  For detection, this
is done with the function \code{\link{io.glm}} which uses an iterative-offset glm as
described by Buckland et al.(1992).  The detection parameters are for the
logit of detection probability.  For pod size parameters, the intial value 0.4 is
used for all but this could be changed to use \code{\link{fit.fS}}.

In the data file, 4 fields must be named specifically.  One is \code{station} which
implies the observer role with \code{P} being primary and \code{S} being secondary.
The records must be in pairs with the record for secondary following the record for
primary.  Another specific field is \code{seen} which is 1 if observer saw the pod
and 0 otherwise.  Also, \code{podsize} must be used as the name of the size of the
pod and \code{Start.year} is the survey year value which is the 4-digit year for
December at the start of the survey.

Also, for the detection formula, use the field \code{True} to specify that pod size
affects detection probability.  There will be no field named \code{True} in the
data but it will be created as part of the estimation procedure.

\code{fit.missed.pods} calls the following functions in the package: \code{\link{lnl.missed.pods}}, 
\code{\link{io.glm}}, \code{\link{podsize.computations}}, and \code{\link{fit.pods}}.  I

\code{lnl.missed.pods} calls contains the code for equations 13-16 in Laake et al.(2009). 

}
\value{
\code{lnl.missed.pods}: returns negative log-likelihood.

\code{fit.missed.pods}: returns a list with elements:
   \item{par}{vector of final parameter estimates}
   \item{AIC}{Akaike information criterion value for model}
   \item{model}{list of model output from optim}
}
\author{Jeff Laake}
\seealso{\code{\link{lnl.missed.pods}},\code{\link{io.glm}}}
