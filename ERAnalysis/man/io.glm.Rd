\name{io.glm}
\Rdversion{1.1}
\alias{io.glm}
\title{Fits detection probability function for double-observer data using iterative
offset approach}
\description{
Code taken from the \code{mrds} package (mark-recapture distance sampling).  It
uses the iterative offset approach described by Buckland et al.(1993) which
allows standard logistic regression to be used for capture history data of
the form 01,10,11 which excludes the 00 which are those missed by both observers.
The conditioning on the set observed is handled through the offset which is
adjusted iteratively based on the estimated detection probabilities until they
converge.
}
\usage{
io.glm(datavec, fitformula, eps = 1e-05, iterlimit = 500)
}
\arguments{
  \item{datavec}{dataframe with the field \code{station} as the
                  observer role}
  \item{fitformula}{formula for the logistic regression}
  \item{eps}{convergence criterion}
  \item{iterlimit}{ maximum number of iterations allowed}
}
\value{
  \item{glmobj}{glm model output}
  \item{comp2 }{final offsetvalues from iterative fit}
}
\references{
BUCKLAND, S. T., J. M. BREIWICK, K. L. CATTANACH and J. L. LAAKE. 1993. Estimated population size of the California gray whale. Marine Mammal Science 9: 235-249.

BUCKLAND, S. T., D. R. ANDERSON, K. P. BURNHAM, J. L. LAAKE, D. L. BORCHERS and L. THOMAS. 2004. Advanced distance sampling: estimating abundance of biological populations. Pages 416. Oxford University Press, Oxford, UK; New York. See Chapter 11.
}
\author{Jeff Laake}
