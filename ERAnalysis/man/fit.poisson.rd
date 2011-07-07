\name{fit.poisson.re}
\Rdversion{1.1}
\alias{fit.poisson.re}
\alias{fit.poisson}
\alias{repois.lnl}
\alias{pois.lnl}
\alias{pois.d}
\title{Fitting code for truncated Poisson distribution with pod size calibration data}
\description{Fits fixed and mixed effects models with a single random effect to pod size calibration data.}
\usage{
fit.poisson(x,formula,initial,maxit=1000,method="BFGS",hessian=FALSE,nmax=20)
fit.poisson.re(x,re.factor,formula,sigma.formula=~1,initial,maxit=1000,method="BFGS",hessian=FALSE,z=5,nmax=20)
repois.lnl(par,mm,ss,observed.size,z=5,nmax=20)
pois.lnl(par,xmat,observed.size,nmax=20)
pois.d(lambda,x,nmax)
}
\arguments{
  \item{x}{pod size calibration dataframe except in \code{pois.d} where it is vector of observed pod sizes}
  \item{formula}{formula for fixed effects portion of model for lambda in Poisson}
  \item{xmat}{model design matrix for lambda constructed from formula and data}
  \item{re.factor}{list of one or more vectors of factor variables to split x into a list of dataframes. Length of each factor variable vector must
                  match number of rows in x}
  \item{sigma.formula}{formula for sigma of normal error}
  \item{observed.size}{vector of observed sizes for \code{pois.lnl}; or \code{repois.lnl} a list of observed size vectors where list elements are split by random effect factor levels}
  \item{initial}{required vector of initial parameter values}
  \item{maxit}{maximum number of iterations (passed to optim)}
  \item{method}{method used by optim for optimization}
  \item{hessian}{passed to optim, if TRUE returns hessian for v-c matrix}
  \item{mm}{list of model matrices for lambda; list elements are split by random effect factor levels}
  \item{ss}{list of model matrices for sigma; list elements are split by random effect factor levels}
  \item{par}{parameter values for model}
  \item{z}{defines integration bounds for normal (default is 5 for integration range of -5*sigma,5*sigma)}
  \item{lambda}{vector of lambda values for observed sizes in x for \code{pois.d}}
  \item{nmax}{maximum pod size}
}
\details{

\code{pois.d} computes the Poisson probabilities for a vector of lambdas matched with observed values (x).
 Because these are pod size estimates (1,2,...) the observed x is shifted (x-1) such that it
 matches a Poisson with values 0,1,2,...  Also, to avoid dealing with an infinite number of
 of possible values the maximum pod size is set (nmax) and the distribution is renormalized
 over the range 0,1,2,...,(nmax-1). See eqn 9 in Laake et al. (2009) for Poisson construction.

\code{pois.lnl} computes negative log-likelihood for a fixed effects Poisson model for
 the pod size calibration data. See eq 8 & 9 in Laake et al (2009).

\code{repois.lnl} computes negative log-likelihood for a mixed effects Poisson model with a single
random effect for the pod size calibration data. A normal 0,sigma error distribution is used
for the random effect.  It is used for the log(lambda) where lambda is the Poisson intensity. See 
equation 10 in Laake et al. (2009) except modified for Poisson.

\code{fit.poisson} fits a Poisson fixed-effects model to the pod size calibration data. Uses \code{pois.lnl}.

\code{fit.poisson.re} fits a Poisson mixed-effects model with a single random effect to the pod size calibration data.
Uses \code{repois.lnl}.

}
\value{
\code{pois.d} returns vector of probabilities 

\code{pois.lnl} and \code{repois.lnl} return negative log-likelihood value.

\code{fit.poisson} and \code{fit.poisson.re} return results list from \code{optim} for
fitted model.
}
\author{Jeff Laake}
\seealso{\code{\link{fit.gamma}}}
\examples{
#
# Fits set of models for pod size calibration data for Poisson distribution in Laake et al.(2009)
# with model selection results shown in Table 4.
#
data(PodsizeCalibrationData)
PodsizeCalibrationData$size=factor(cut(PodsizeCalibrationData$True,c(1,2,3,4,21),right=FALSE))
levels(PodsizeCalibrationData$size)=c("1","2","3","4+")
PodsizeCalibrationData$plussize=as.numeric(PodsizeCalibrationData$True>=4)
poisson.observer=fit.poisson.re(x=PodsizeCalibrationData,re.factor=PodsizeCalibrationData$Observer,formula=~size+True:plussize,
                                initial=c(-1.5639437 ,-1.5340573 , 1.1582312 , 1.6418947 , 1.1796440,  0.2231416))

poisson.pod=fit.poisson.re(x=PodsizeCalibrationData,re.factor=factor(paste(PodsizeCalibrationData$Year,PodsizeCalibrationData$ID)),formula=~size+True:plussize,
                                initial=c(-0.8054231, -1.6909508,  1.3141876 , 1.8451641 , 1.6701502 , 0.1498825))

poisson.pod.size=fit.poisson.re(x=PodsizeCalibrationData,re.factor=list(factor(paste(PodsizeCalibrationData$Year,PodsizeCalibrationData$ID)),PodsizeCalibrationData$size),formula=~size+True:plussize,
                       initial=c(-0.8054231,0,0,0, -1.6909508,  1.3141876 , 1.8451641 , 1.6701502 , 0.1498825),
                       sigma.formula=~size)

poisson.year=fit.poisson.re(x=PodsizeCalibrationData,re.factor=PodsizeCalibrationData$Year,formula=~size+True:plussize,
                                initial=c(-1.5206800 ,-1.4598787,  1.1269417 , 1.6181561,  1.2654947 , 0.2070926))

poisson.dot=fit.poisson(PodsizeCalibrationData,formula=~size+True:plussize,
              initial=c(-1.690892, 1.2608829, 1.625358, 1.6072374, .345172),maxit=2000)

poisson.year.fixed=fit.poisson(PodsizeCalibrationData,formula=~Year*(size+True:plussize),
              initial=c(-1.690892, 0.2608829, 0.0625358, 0.6072374, 1.345172, 1.449742, 0.7439905 ,0.4038016, -0.3388784 ,
              -0.2160762, -0.6082892, 0.1826372, 0.5286403, 0.3515043, -0.57903, 0.05635473, -0.6360614, 0.1605690, -0.1755681 ,-0.0709657),
              maxit=3000)
AICc=function(model,n)
{
  K=length(model$par)
  return(2*model$val+2*K*(n/(n-K-1)))
}

AICc(poisson.dot,196)
AICc(poisson.year.fixed,196)
AICc(poisson.pod,196)
AICc(poisson.observer,196)
AICc(poisson.year,196)

}