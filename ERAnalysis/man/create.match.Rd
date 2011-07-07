\name{create.match}
\Rdversion{1.1}
\alias{create.match}
\alias{matchdf}
\alias{linklist}
\alias{linkdf}
\alias{CreateObserverExperience}
\title{Linking and matching algorithms for gray whale double observer data}
\description{
Only \code{create.match} is intended to be invoked by the user.  It in turn calls the
other functions to optionally link (\code{lcut >0}) within an observer record and then match across observer records.

\code{linklist} uses function \code{linkdf} to identify potential links within each observer record (P: primary and S:secondary) and creates 
a revised dataframe using the averages of t241 and distance and the summed group size.  

\code{matchdf} finds matches of the same gray whale pods seen by the observers at the
P and S stations. 
}
\usage{
match.data=create.match(lcut=0.2,mcut=1,twt=.18,dwt=3.95,pwt=0.05,crittype="ellipse",fdir="c:/gw")
linklist(mspdf,lcut=-0.05,twt=.18,dwt=3.95,crittype="ellipse")
linkdf(mspdf,lcut=-0.05,twt=.18,dwt=3.95,crittype="ellipse")
matchdf(mspdf,mcut=1,twt=.18,dwt=3.95,pwt=.05,crittype="ellipse")
CreateObserverExperience()
}
\arguments{
  \item{mspdf}{dataframe of sighting records with specific field station to distinguish between primary and secondary stations}
  \item{lcut}{linking parameter cutoff}
  \item{mcut}{matching parameter cutoff}
  \item{twt}{time weight for crossing time difference in minutes}
  \item{dwt}{distance weight for ratio of difference to longer distance of shore at 241}
  \item{pwt}{weight for difference in pod size}
  \item{crittype}{criterion for linking: "box" criterion is maximum of distances; "diamond" criterion is sum of distances and "ellipse" 
       criterion is the euclidean distance}
  \item{fdir}{directory containing ERAbund output files}
}
\details{
Finds matches (same pod seen by both observers) for double observer data based 
on arguments for the weights, criterion and cutoffs. If \code{lcut>0} pods seen by an observer that are in close proximity
are linked and merged into combined pods.  This is done for both observers and then the linked pods
are compared across observers for matches.  See Appendix of Laake et al. for a more complete description.
It uses \code{\link{linkdf}} and \code{\link{matchdf}}.

ERAbund: Uses the diamond spatial criterion with user specified parameters,
         default is the criteria used by Breiwick et all 2006

 0.04   !Time weight factor for linking & matching

 1.75   !Distance weight factor for linking & matching

 0.00   !Pod size weight factor for linking

 0.05   !Pod size weight factor for matching

-0.05   !Default maximum sum score for linking (merging pods) - which means no linking

 1.0    !Default maximum sum score for matching
}
 
\value{
From \code{create.match} returns a dataframe with 2 records (one for north and one for south station) for
each uniquely determined gray whale pod.  The field \code{seen} is 1 if it was seen by the observer at the
station and 0 if not.  See \code{\link{Match}} for a description of fields.

From \code{CreateObserverExperience}
\item{ObserverExp}{dataframe of cummulative hours of experience for each observer for each survey date from 1967-2007}

}
\author{Rod Hobbs, Jeff Laake}



