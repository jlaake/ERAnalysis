\name{compute.series}
\Rdversion{1.1}
\alias{compute.series}
\title{Compute gray whale abundance for 23 year series}
\description{
For a set of detection models, fit best model for each year from 1987-2006, compute abundance for each
year, then form conventional ratio estimator with best and naive abundance estimates and use ratio to compute abundance
estimates for 1967-1984. 
}
\usage{
compute.series(models,naive.abundance,sightings=NULL,effort=NULL,
          gsS=NULL,best=TRUE,Match=NULL,cutoff=4,final.time=c(rep(90,20),100,100,90),
          lower.time=rep(0,23),lcut=0.2,mcut=1.0,twt=0.18,dwt=3.95,pwt=0.05,crittype="ellipse",
          DistBreaks=c(0,1,2,3,4,20),Use=TRUE,hessian=FALSE,debug=FALSE,TruePS=TRUE,
          recent.years=c(1987,1992,1993,1995,1997,2000,2001,2006),fn=1.0817)
}
\arguments{
  \item{models}{vector of character strings specifying set of detection models to examine}
  \item{naive.abundance}{set of naive abundance models and estimates for all 23 years; the naive estimates do not
                           include corrections for detection and pod size -- see example}
  \item{sightings}{primary observer sightings to be used in the analysis; if NULL uses PrimarySightings}
  \item{effort}{primary observer effort to be used in the analysis; if NULL uses PrimaryEffort}
  \item{gsS}{pod size calibration matrix; if NULL it retrieves stored dataframe}
  \item{best}{if TRUE, only computes abundance for best detection model (lowest AIC)}
  \item{Match}{dataframe of matched double observer records; if NULL, the dataframe is constructed using specified parameters}
  \item{cutoff}{cutoff for delta AIC to select set of detection models}
  \item{final.time}{vector of final times to use for anchoring}  
  \item{lower.time}{vector of beginning times to use for anchoring}
  \item{lcut}{linking parameter cutoff}
  \item{mcut}{matching parameter cutoff}
  \item{twt}{time weight for crossing time difference in minutes}
  \item{dwt}{distance weight for ratio of difference to longer distance of shore at 241}
  \item{pwt}{weight for difference in pod size}
  \item{crittype}{criterion for linking: "box" criterion is maximum of distances; "diamond" criterion is sum of distances and "ellipse" 
       criterion is the euclidean distance}
  \item{DistBreaks}{breaks in offshore distance to create distance factor variable for detection}
  \item{Use}{if TRUE, only use watches that meet the criterion used to select pre-1987 data of beaufort<4 and vis<=4}
  \item{hessian}{if TRUE, compute and return the hessian for the detection and pod size parameters}
  \item{debug}{if TRUE, print intermediate results when fitting detection model and pod size distribution}
  \item{TruePS}{if TRUE, use of True is forced in models for pod size; this takes account unknown true pod size in detection function}
  \item{recent.years}{vector of recent survey Start.year values}
  \item{fn}{nighttime multiplicative correction factor to be applied to series}
}
\details{
The following steps are followed in this function:

1) The sightings are merged with their matching effort record and if \code{Use==TRUE}, then only those
sightings and effort with \code{Use==TRUE} are included in the analysis.  This eliminates any effort
or sighting records that occurred during a watch in which beaufort or vis exceeded 4 (based on the
primary observer's record).  This is done to be consistent with filtering of records for surveys
prior to 1987.

2) Next if lcut > 0, the remaining primary sightings are matched against each other to see if
any should be merged based on their proximity in terms of distance and passage time (t241).  "Proximity" is
determined by the value of lcut and the weighting factors twt (for time) and dwt (for distance). If they 
are linked, for the merged sightings,the distances are averaged and the observed pod sizes are summed.

3) Next for all years from 1987-2006, all sightings when there were 2 observers are matched across observer with the
function \code{\link{create.match}} with the same lcut parameters and the same values used for twt and dwt
and with the additional pwt (for pod size) used in the matching criterion. These matched data which have 2 records for each pod seen (by either observer) with one for the
primary observer and the other for the secondary observer (station="P" or "S").  The field "seen" is either
0 or 1 depending on whether it was missed or seen by that observer.  There are also other fields such
as distance, vis, beaufort, pod size etc.

4) If \code{Use==TRUE}, the matched data are merged with the effort records and any with the field \code{Use==FALSE}
in the effort are excluded.  The value of Use is determined based on the primary observer's assessment
of beaufort and vis as it is the primary observer's data that are used to estimate abundance. In some cases,
the matched data will not have an effort records because the effort has already been excluded because it has
already been filtered for vis and beaufort.  These are also excluded.

5) Models fitted with \code{\link{fit.missed.pods}} take a long time to run and need starting values.
To provide starting values and to select the best model we use an approximate method to analyze the
double observer data which assumes no error in pod size. For each of the 8 surveys from 1987-2006 (\code{recent.years}) that had double observer data, the specified set of 
detection models are fitted using the iterative offset glm (\code{\link{io.glm}})and the best set of models is 
selected for each year.  For this initial model fitting and selection, the observed pod 
size is used instead of the unknown true pod size.  The set of models retained
is determined by the value of \code{cutoff}.  Only models with a delta AIC less than cutoff is
retained.  Once that set is selected, it is further refined by removing any models that have
positive estimates for beaufort or visibility which is counter to what they should be.  A positive
value for either of these factors most likely indicates a missing covariate for the model.

6) For each of the 8 recent surveys, the detection parameters and pod size parameters are fitted 
with \code{\link{fit.missed.pods}}.  Those estimated parameters and data are fed to \code{\link{estimate.abundance}}
to estimate abundance using a Horvitz-Thompson estimator for predicted number of
whales from the observations.  With this series of abundance predictions, a gam smooth over time is
fitted and used to predicted the total abundance of whales that migrated past that season. There are
2 potential calls for abundance estimation to handle estimation for the previous approach which 
used observed pod size in the detection function and corrected pod size with the additive correction factor.
If \code{TruePS==TRUE} or there is no \code{corrected.podsize} field in the data,
then the newer approach is implemented with the true pod size parameters determined with the fitting of the 
detection parameters.  If there is a \code{corrected.podsize} field in the data and \code{TruePS==FALSE},
then the abundance is estimated using the corrected pod size values.

7) Once abundance is estimated for each of the 8 recent surveys, a ratio estimator is constructed with the abundance
estimate and the naive abundance estimate for each of the 8 years (argument \code{naive.abundance}. 
The naive abundance estimates for the 15 surveys prior to 1987 are multiplied by the ratio to construct an estimate with "correction" for
missed pods and pod size bias.

8) The final set of abundance estimates (Nhat) include an adjustment (fn) for a differential rate of migration
during the latter half of the migration. 

This function calls the following functions in the package: \code{\link{linkdf}}, \code{\link{estimate.abundance}},
 \code{\link{io.glm}}, \code{\link{create.match}}, \code{\link{fit.missed.pods}}.  It contains the code for
equations 23,24,and 29 in Laake et al.(2009).
}
\value{
\item{summary.df}{dataframe of summary statistics for last 8 surveys}
\item{Match}{dataframe of matched data for double observers}
\item{formulae}{list of selected detection model formula for each year}
\item{detection.models}{list (for each of 8 years) of list of fitted detection/pod size models}
\item{abundance.models}{list (for each of 8 years) of list of fitted abundance models and estimates (without fn)}
\item{ratio}{abundance/naive abundance ratio estimate}
\item{Nhat}{vector of abundance estimates for all 23 years which includes nighttime correction factor}
}
\author{Jeff Laake}



