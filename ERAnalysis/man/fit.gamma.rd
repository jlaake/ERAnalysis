\name{fit.gamma}
\Rdversion{1.1}
\alias{fit.gamma.re}
\alias{fit.gamma}
\alias{regam.lnl}
\alias{gam.lnl}
\alias{gam.d}
\title{Fitting code for discretized gamma distribution with pod size calibration data}
\description{Fits fixed and mixed effects gamma models with a single random effect to pod size calibration data.}
\usage{
fit.gamma.re(x,re.factor,sformula,rformula,sigma.formula=~1,initial,maxit=1000,method="BFGS",hessian=FALSE,nmax=20,z=5)
regam.lnl(par,mm,rr,ss,observed.size,nmax=20,z=5)
fit.gamma(x,sformula,rformula,initial,maxit=1000,method="BFGS",hessian=FALSE,nmax=20)
gam.lnl(par,smat,rmat,observed.size,nmax=20)
gam.d(shape,rate,x,nmax)
}
\arguments{
  \item{x}{pod size calibration dataframe except in \code{gam.d} where it is vector of observed pod sizes}
  \item{re.factor}{list of one or more vectors of factor variables to split x into a list of dataframes. Length of each factor variable vector must
                  match number of rows in x}
  \item{sformula}{formula for fixed effects portion of model for shape in gamma}
  \item{rformula}{formula for fixed effects portion of model for rate in gamma}
  \item{sigma.formula}{formula for sigma of normal error}
  \item{smat}{model design matrix for shape constructed from \code{sformula} and data}
  \item{rmat}{model design matrix for rate constructed from \code{rformula} and data}
  \item{observed.size}{vector of observed sizes for \code{gam.lnl}; or \code{regam.lnl} a list of observed size vectors where list elements are split by random effect factor levels}
  \item{initial}{required vector of initial parameter values}
  \item{maxit}{maximum number of iterations (passed to optim)}
  \item{method}{method used by optim for optimization}
  \item{hessian}{passed to optim, if TRUE returns hessian for v-c matrix}
  \item{mm}{list of model matrices for rate; list elements are split by random effect factor levels}
  \item{rr}{list of model matrices for shape; list elements are split by random effect factor levels}
  \item{ss}{list of model matrices for sigma; list elements are split by random effect factor levels}
  \item{par}{parameter values for model}
  \item{z}{defines integration bounds for normal (default is 5 for integration range of -5*sigma,5*sigma)}
  \item{shape}{vector of shape values for observed sizes in x for \code{gam.d}}
  \item{rate}{vector of rate values for observed sizes in x for \code{gam.d}}
  \item{nmax}{maximum pod size}
}
\details{

\code{gam.d} computes the gamma probabilities for a vector of lambdas matched with observed values (x).
 Because these are pod size estimates (1,2,...) the observed x is shifted (x-1) such that it
 matches a gamma with values 0,1,2,...  Also, to avoid dealing with an infinite number of
 of possible values the maximum pod size is set (nmax) and the distribution is renormalized
 over the range 0,1,2,...,(nmax-1). See eqn 7 & 10 in Laake et al. (2009) for gamma construction.

\code{gam.lnl} computes negative log-likelihood for a fixed effects gamma model for
 the pod size calibration data. See eq 8 & 9 in Laake et al but modified for gamma model.

\code{regam.lnl} computes negative log-likelihood for a mixed effects gamma model with a single
random effect for the pod size calibration data. A normal 0,sigma error distribution is used
for the random effect.  It is used for the log(lambda) where lambda is the gamma intensity. See 
equation 10 in Laake et al. (2009).

\code{fit.gamma} fits a gamma fixed-effects model to the pod size calibration data. Uses \code{gam.lnl}.

\code{fit.gamma.re} fits a gamma mixed-effects model with a single random effect to the pod size calibration data.
Uses \code{regam.lnl}.

}
\value{
\code{gam.d} returns vector of probabilities

\code{gam.lnl} and \code{regam.lnl} return negative log-likelihood value.

\code{fit.gamma} and \code{fit.gamma.re} return results list from \code{optim} for
fitted model.
}
\author{Jeff Laake}
\seealso{\code{\link{fit.poisson}}}
\examples{
#
# Fits set of models for pod size calibration data for gamma distribution in Laake et al.(2009)
# with model selection results shown in Table 4.
#
data(PodsizeCalibrationData)
PodsizeCalibrationData$size=factor(cut(PodsizeCalibrationData$True,c(1,2,3,4,21),right=FALSE))
levels(PodsizeCalibrationData$size)=c("1","2","3","4+")
PodsizeCalibrationData$plussize=as.numeric(PodsizeCalibrationData$True>=4)

gamma.dot=fit.gamma(PodsizeCalibrationData,sformula=~size+True:plussize,rformula=~size,
              initial=c(-.128,1,1,.9,-.95,.44,.5,.2,-.3),maxit=2000)

gamma.year.fixed=fit.gamma(PodsizeCalibrationData,sformula=~Year*(size+True:plussize),rformula=~Year*size,
              initial=rep(0,36),maxit=2000)


gamma.observer=fit.gamma.re(PodsizeCalibrationData,re.factor=list(PodsizeCalibrationData$Observer,PodsizeCalibrationData$size),
                                rformula=~size+True:plussize,sformula=~size,sigma.formula=~1,
                                initial= c(-1.590104,0.8107313,0.5089849,-0.03983425,0.2423928,-0.2034048,0.2830661,1.189406,0.9777272,0.7808261),maxit=2000)


gamma.pod=fit.gamma.re(PodsizeCalibrationData,re.factor=list(factor(paste(PodsizeCalibrationData$Year,PodsizeCalibrationData$ID)),PodsizeCalibrationData$size),
                                rformula=~size+True:plussize,sformula=~size,sigma.formula=~1,
                                initial= c(-0.9214691,1.184905,0.4498295,0.003455527,0.3922958,-0.1893637,0.683027,1.071157,1.092008,0.8718244),maxit=2000,
                                method="Nelder-Mead")

gamma.year=fit.gamma.re(PodsizeCalibrationData,re.factor=list(PodsizeCalibrationData$Year,PodsizeCalibrationData$size),
                                rformula=~size+True:plussize,sformula=~size,sigma.formula=~1,
                                initial= c(-0.9214691,1.184905,0.4498295,0.003455527,0.3922958,-0.1893637,0.683027,1.071157,1.092008,0.8718244),maxit=2000,
                                method="Nelder-Mead")

gamma.pod=fit.gamma.re(PodsizeCalibrationData,re.factor=list(factor(paste(PodsizeCalibrationData$Year,PodsizeCalibrationData$ID)),PodsizeCalibrationData$size),
                                rformula=~size+True:plussize,sformula=~size,sigma.formula=~1,
                                initial= gamma.pod$par,maxit=2000,hessian=TRUE,method="Nelder-Mead")

AICc=function(model,n)
{
  K=length(model$par)
  return(2*model$val+2*K*(n/(n-K-1)))
}

AICc(gamma.dot,196)
AICc(gamma.year.fixed,196)
AICc(gamma.pod,196)
AICc(gamma.observer,196)
AICc(gamma.year,196)
}