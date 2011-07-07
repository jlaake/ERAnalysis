\name{gamma.pod}
\alias{gamma.pod}
\docType{data}
\title{Fitted Gamma Pod Size Calibration Model with Random Pod Effect}
\description{
Results list from fitted model used to create \code{\link{gsS}}.  It contains the 
parameters and hessian for variance estimation.}
\usage{
data(gamma.pod)
}
\format{             
  \describe{
    \item{\code{gamma.pod}}{A list of results from call to optim for model fitted to 
                            pod size calibration data}
    }
}
\details{
The fitted model and resulting hessian are a bit dicey and can hinge on the starting values. To
document the fitted model used to obtain the results it has been saved here.  Hessians from optim
are not always numerically stable. The hessian function in the numderiv package is better but when
I tried that function it wanted more memory than was avaialable (>2GB).  For the future, it may be
worthwhile investigating a bootstrap technique to estimate the variance-covariance matrix for the parameters
of the model.
}
\keyword{datasets}
