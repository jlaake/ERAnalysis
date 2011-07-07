\name{compute.series.var}
\Rdversion{1.2}
\alias{compute.series.var}
\title{Computes variance-covariance matrix for gray whale abundance for 23 year series}
\description{
Computes the variance-covariance matrix for the estimated abundances for the 23 year time series.
There are three components for the variance: 1) variation in predicted abundances from estimation of detection and pod size parameters,
2) variation in predicted abundances from estimation of pod size calibration matrix,and 
3) variation in predicted abundances from fitted migration rate gam and residual variation of migration rate 
}
\usage{
compute.series.var(results,naive.abundance,ps.results,Use=TRUE,lcut=0.2,twt=0.18,dwt=3.95,
                   crittype="ellipse",final.time=c(rep(90,20),100,100,90), 
                   lower.time=rep(0,23), DistBreaks=c(0,1,2,3,4,20),
                   recent.years=c(1987,1992,1993,1995,1997,2000,2001,2006),  
                   fn=1.0817,se.fn=0.0338, debug=FALSE, delta=c(0.001,0.01))
}
\arguments{
  \item{results}{list of results from \code{\link{compute.series}}}
  \item{naive.abundance}{set of naive abundance models and estimates for all 23 years; the naive estimates do not
                           include corrections for detection and pod size -- see example}
  \item{ps.results}{list with the first element being the fitted \code{gamma.pod} model and the second named vc 
                     containing the var-cov matrix for the parameters in the model}
  \item{Use}{if TRUE, only use watches that meet the criterion used to select pre-1987 data of beaufort<4 and vis<=4}
  \item{lcut}{linking parameter cutoff}
  \item{twt}{time weight for crossing time difference in minutes}
  \item{dwt}{distance weight for ratio of difference to longer distance of shore at 241}
  \item{crittype}{criterion for linking: "box" criterion is maximum of distances; "diamond" criterion is sum of distances and "ellipse" 
       criterion is the euclidean distance}
  \item{final.time}{vector of final times to use for anchoring}  
  \item{lower.time}{vector of beginning times to use for anchoring}
  \item{DistBreaks}{breaks in offshore distance to create distance factor variable for detection}
  \item{recent.years}{vector of recent survey Start.year values}
  \item{fn}{nighttime multiplicative correction factor to be applied to series}
  \item{se.fn}{standard error for nighttime multiplicative correction factor to be applied to series}
  \item{debug}{if TRUE will output some intermediate results for derivative calculation for var2 component}  
  \item{delta}{vector of two proportional changes for numerical derivatives for variance component 1 and 2}  
}
\details{
For variance components 1 and 2, a numerical delta method approximation for the variances and covariances
are computed.  This involves computing the abundances at +/- delta value away from the maximum likelihood
estimate for each parameter to estimate a first derivative for each of the abundances for each
parameter.  The variance-covariance matrix of the parameters is then pre- and post-multiplied by the
numerical derivative matrix to obtain the variance-covariance matrix of the abundances for the set of
estimated model parameters.  

The first variance component is from the estimation of the parameters of the
gamma distribution for true pod size and for the parameters of the detection model. The numerical
derivative matrix for this component is easily calculated by adjusting the model parameter (+/- delta) and
calling \code{\link{estimate.abundance}} to refit the migration gam and compute total abundance.  The
pod size and detection parameters were fitted separately for each of the recent survey years (1987-2006); thus,
for this component there is no covariance between the most recent 8 surveys.  However,
the average ratio of the estimated abundance to naive abundance for the 8 years is used to "correct" the
naive abundance estimates prior to 1987 which induces a correlation between the abundance estimates prior to 1987
and those for 1987 and beyond.  

The second component of variance is due to estimation of the parameters of the pod size calibration matrix
which is used in both fitting the detection/pod size parameters which in turn affects the estimated
number of whales passing and the total abundance.  Thus, to compute the numerical derivatives for
this component it is necesssary to re-fit the detection/pod size model for each of the 8 recent years for
each adjustment of the 8 pod size calibration parameters.  This can take a long time to compute -- something on
the order of 6 hours and a relatively fast machine.  Because a common pod size calibration matrix is
used for all of the abundance estimates, this creates a covariance between each pairing of the 23 estimates.
Beginning in v1.5, the proportional change was increased to 0.01 from 0.001 for derivative calculation, so
the change in gsS was sufficiently large so as not to be diluted by error level in optimization.

The final component of variance is much easier to compute. See equations 25-28 in Laake et al. (2009).
This component captures the variation due to fitting the gam to the migration rate which is used to compute 
total abundance from predicted abundance in each survey period.  This component is computed conditional on 
the estimated values of the detection/pod size parameters and the parameters of the pod size calibration matrix. The variances and 
covariances are computed from the fitted gam and the residual variation of the observed rate around the
gam predicted rate.  Variances for the older estimates (<1987) and their covariance with the more recent estimates (>=1987) due to
ratio correction factor are included as well using the variances from the fitted gam for the naive 
abundance estimates.

The total variance-covariance matrix is the sum of the 3 component variances.  The total v-c matrix
is then adjusted for the nighttime correction factor.  The correction factor is fn=1.0875 (Rugh et al. 2005)
with an std error of 0.0363.  The covariances are multiplied by fn^2 and the variances are adjusted as
var(N*fn)=(N*fn)^2*(cv^2(fn)+cv^2(N)). See equations 29-31 in Laake et al. (2009). 

The standard errors of the abundance estimates is the square root of the diagonal elements.  The
correlation matrix for the abundance estimates is also returned.

This function calls the following functions in the package: \code{\link{linkdf}}, \code{\link{estimate.abundance}},
 \code{\link{create.gsS}}, \code{\link{fit.missed.pods}}.  It contains the code for
equations 23 to 24 and 27 to 31 in Laake et al.(2009).  Results for Eq 25&26 are obtained from \code{\link{fit.migration.gam}} which
is called by \code{estimate.abundance}.
}
\value{
\item{vc}{total variance-covariance matrix of 23 abundance estimates}
\item{var1}{variance-covariance matrix of 23 abundance estimates due to fitted detection/pod size parameters}
\item{var2}{variance-covariance matrix of 23 abundance estimates due to fitted pod size calibration parameters}
\item{var3}{variance-covariance matrix of 23 abundance estimates due to fitted gam for migration rate and residual variation in rate}
\item{se}{standard errors of 23 abundance estimates \code{sqrt(diag(vc))}}
\item{cormat}{correlation matrix for 23 abundance estimates}
}
\author{Jeff Laake}



