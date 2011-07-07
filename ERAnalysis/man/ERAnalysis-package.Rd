\name{ERAnalysis-package}
\Rdversion{1.4}
\alias{ERAnalysis-package}
\alias{ERAnalysis}
\docType{package}
\title{
Analysis code and data for gray whale southbound migration survey
}
\description{
Functions for reading/merging/processing data, fitting pod size calibration data, analysis of pods missed during survey, and estimation of
abundance using the various components.
}
\details{
\tabular{ll}{
Package: \tab ERAnalysis\cr
Type: \tab Package\cr
Version: \tab 1.4\cr
Date: \tab 2009-11-05\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}

For abundance estimation, the essential dataframes are: 1)\code{\link{PrimaryEffort}} and 1)\code{\link{PrimarySightings}} for computation of
the encounter rate/migration curve, 2) the matched observer data created by \code{\link{create.match}} which is used to compute a detection probability function for the primary observer from
the capture-recapture data created with the secondary independent observer, and 3) 1)\code{\link{gsS}} the pod size calibration matrix that is computed 
from the pod size calibration data (\code{\link{PodsizeCalibrationData}}, \code{\link{PodsizeCalibrationTable}}) and the code described in 
\code{\link{create.podsize.calibration.matrix}}. The dataframes are contained in the package and can be obtained 
by using the \code{data} function.  For example, \code{data(PrimaryEffort)} will get that dataframe from the
package data directory and put it in the R workspace so it can be used.  The entire list of dataframes available
in the package are listed at the bottom of this help file with a listing and description of the functions for 
analysis. More detail is given in \code{\link{Maintenance}} on how the dataframes were constructed from the various sources of data collected from
the surveys spanning 40 years.  The maintenance information is primarily needed if there are any future edits to the data in the database
to incorporate the modified data into the package. 


A brief description of each analysis function is given below with a link to a more complete help file:

Pod size calibration analysis:

1) \code{\link{create.podsize.calibration.matrix}}: creates \code{\link{gsS}} the pod size calibration matrix from the 
the fitted gamma model using a random pod effect (\code{\link{create.gsS}}) and also computes the additive correction factor vectors \code{\link{add.cf.all}}, 
\code{\link{add.cf.reilly}}, and \code{\link{add.cf.laake}} and
saves all 4 data sets in the package data directory.
   
2) \code{\link{fit.gamma}}: fits a truncated discretized gamma distribution to
the observed pod sizes with fixed effects model.   

3) \code{\link{fit.poisson}}: fits a truncated Poisson distribution to
the observed pod sizes with a fixed effects model.     

4) \code{\link{fit.gamma.re}}: fits a truncated discretized gamma distribution to
the observed pod sizes with mixed effects model with a single random effect.   

5) \code{\link{fit.poisson.re}}: fits a truncated Poisson distribution to
the observed pod sizes with  mixed effects model with a single random effect.     
  
6) \code{\link{fit.fS}}: fits a gamma distribution for true pod size from the input
pod size calibration matrix \code{\link{gsS}} and the distribution of observed pod sizes. 

7) \code{\link{reilly.cf}}: create 'corrected' pod sizes using the approach first suggested by Reilly.

8) \code{\link{expected.podsize}}: computes the true pod size distribution \code{fS} and the expected pod size can be computed from the gamma distribution
specified by the parameter argument are computed as well as the conditional distribution of the true size \code{S} given an observed size \code{s} 
(\code{fSs}) and conditional expectation \code{ESs} based on the pod size calibration matrix \code{gsS}

 
Migration curve analysis:

1) \code{\link{fit.migration.gam}}: Fits a GAM to the migration data over time using either number of pods or number
of whales.  It can optionally anchor the data to include 0 observations at times \code{lower.time} and \code{final.time}.
This function is called from \code{\link{estimate.abundance}} but can be called on its own with the proper
arguments. It is used for estimation of naive abundance and for abundance corrected for missed pods and
pod size error. 

2) \code{\link{compute.sampling.multiplier}}:Integrates the area under the fitted GAM migration curve and the area under
the curve contained in the intervals that were sampled.  The ratio of the
total area to the area within the intervals can be used as a multiplier to expand
the abundance for the entire migration; however, \code{fit.migration.gam} uses a sum of 
daily predictions which is easier for variance computation.

Double-observer analysis:

1) \code{\link{extract.match.data}}: extracts data for matching during double-observer
surveys (experiment=2).  Creates a list of dataframes where each element is the data for
a survey date that contained one or more watches of experiment=2.

2) \code{\link{create.match}}: Finds matches (same pod seen by both observers) for double observer data based 
on arguments for the weights, criterion and cutoffs. If \code{lcut>0} pods seen by an observer that are in close proximity
are linked and merged into combined pods.  This is done for both observers and then the linked pods
are compared across observers for matches.  See Appendix of Laake et al. for a more complete description.
It uses \code{\link{linkdf}} and \code{\link{matchdf}}.

3) \code{\link{io.glm}}:Code taken from the \code{mrds} package (mark-recapture distance sampling).  It
uses the iterative offset approach described by Buckland et al. (1993) which
allows standard logistic regression to be used for capture history data of
the form 01,10,11 which excludes the 00 which are those missed by both observers.
The conditioning on the set observed is handled through the offset which is
adjusted iteratively based on the estimated detection probabilities until they
converge. This approach assumes pod size is known without error and is used to produce starting 
values for and to select the model used in \code{\link{fit.missed.pods}}.

4) \code{\link{fit.missed.pods}}:Fits true pod size distribution and the detection probability parameters
for a given detection model specified by formula. 

5) \code{\link{lnl.missed.pods}}:  Computes negative log-likelihood for true pod size distribution and
detection probability parameters for double observer survey data. Called by \code{\link{fit.missed.pods}}. 

6) \code{\link{fit.pods}}:Fits true pod size gamma distribution conditional on a set of input detection probability parameters.  
and the pod size calibration matrix \code{\link{gsS}} to the distribution of observed pod sizes for a single observer.

7) \code{\link{fit.pods}} Computes negative log-likelihood at given estimates parameters of gamma distribution 
for true pod size from counts of recorded pod sizes (for a single observer) and an estimated calibration matrix.
Called by \code{\link{fit.pods}}.

8) \code{\link{podsize.computations}}: returns a list of fitted gamma distributions for true pod size using \code{\link{fit.fS}}
which assumes pod size does not affect detectability.  These fitted distributions are used for starting values for
the analysis that does allow for pod size to affect detection (\code{\link{fit.pods}} and \code{\link{fit.missed.pods}}). 
It also optionally computes some plots of observed distribution. 


Abundance estimation:

1) \code{\link{estimate.abundance}}: Uses \code{\link{fit.migration.gam}} to fit a migration curve to a sequence of effort periods
and whale passage rates to create naive abundance estimates or abundance estimates
corrected for missed pods and pod size error.  The latter uses a fitted detection function 
and parameters for pod size distribution to construct estimates of abundance from the observations 
using a Horvitz-Thompson estimation approach (eg 1/p where p varies by observation). From the fitted
migration curve it computes total whale abundance.

2) \code{\link{compute.series}}: For a set of detection models, fit best model for each year from 1987-2006, compute abundance for each
year, then form conventional ratio estimator of best and naive abundance estimates and ratio to compute abundance
estimates for 1967-1985.

3) \code{\link{compute.series.var}}: Computes the variance-covariance matrix for the estimated abundances for the 23 year time series.
There are three components for the variance: 1) variation in predicted abundances from estimation of detection and pod size parameters,
2) variation in predicted abundances from estimation of pod size calibration matrix,and 
3) variation in predicted abundances from fitted migration rate gam and residual variation of migration rate 


The following data are available with the package:

1) \code{\link{EarlyEffort}}: all effort for surveys from 1967-1985

2) \code{\link{EarlySightings}}: all sightings for surveys from 1967-1985

3) \code{\link{ERSurveyData}}: all data records for surveys from 1987-2006

4) \code{\link{Observer}}: observers for surveys from 1967-2006

5) \code{\link{Primary}}: on-effort southbound gray whale sightings from Primary observer for surveys from 1987-2006 

6) \code{\link{PrimaryOff}}: off-effort southbound gray whale sightings from Primary observer for surveys from 1987-2006 

7) \code{\link{PrimarySightings}}: on-effort vis and beaufort filtered southbound gray whale sightings from Primary observer for surveys from 1967-2006 

8) \code{\link{SecondarySightings}}: all southbound gray whale sightings from Secondary observer for surveys from 1987-2006 

9) \code{\link{Match}}: matched and unmatched sightings from double-observer survey periods for surveys from 1987-2006 as 
determined by the ERAbund code. These are not used in the abundance estimation. Instead, the \code{\link{create.match}} and accompanying
code is used to create the matches between sightings of the independent observers.

10) \code{\link{PrimaryEffort}} : vis and beaufort filtered effort records for primary observer for all surveys from 1967-2006; \code{AllPrimaryEffort } uses
10 for beaufort and vis limits which will include all effort.

11) \code{\link{SecondaryEffort}}: vis and beaufort filtered effort records for secondary observer for surveys from 1987-2006; \code{AllSecondaryEffort } uses
10 for beaufort and vis limits which will include all effort.

12) \code{\link{gsS}}: pod size calibration matrix

13) \code{\link{PodsizeCalibrationData}}: pod size calibration dataframe

14) \code{\link{PodsizeCalibrationTable}}: pod size calibration data table

15) \code{\link{add.cf.all}},\code{\link{add.cf.reilly}},\code{\link{add.cf.laake}}: vectors of additive pod size correction factors

16) \code{\link{Height}}: observer eye heights for calculation of distance from reticles

Updates:

Version 1.2:
\itemize{
\item  Added nighttime correction factor to abundance estimation which impacted \code{\link{compute.series}} and \code{\link{compute.series.var}}
\item Added more complete example to compute results ERAnalysis.  Code (commented out) will produce a pdf with the
confidence intervals.  To do so the package now requires the Hmisc package.
\item Changed default twt to 0.18 and dwt to 3.95.
}
Version 1.3:
\itemize{
\item Changed to conventional ratio estimator to correct naive abundance estimates prior to 1987
\item Added code to fit random effects model for pod size calibration data and changed saved gsS to
        use the gamma distribution model with the pod random effect
\item Changed default twt to 0.18 and dwt to 3.95 in matching code.
\item Changed plotting routine to use ggplot2 (Hmisc no longer required) and included previous abundance estimates
}

Version 1.4:
\itemize{
\item Fixed some portions of the example code below.
\item Updated and clarified portions of various help files.
}

Version 1.5:
\itemize{
\item Added an argument \code{debug} for \code{\link{compute.series.var}} to output intermediate
results for computation of derivatives in computation of var2 component.  In addition,
I added a delta vector with defaults of 0.001 and 0.01 for the proportional change in the parameters
for computation of the numerical derivatives of var1 and var2 components respectively. For the hi case, a
stable v-c matrix was achieved with delta[2]=0.02 for var2.  In addition, in fitting the models I now
use the fitted model values as initial values and I tightened up the derivative calculations for the optimization in \code{\link{fit.missed.pods}}.
This helped improved derivative calculations because previously the change in the pod size
calibration parameters was getting swamped out by error in the optimization for the
detection/pod parameter estimation.  By tightening up the optimization and increasing the
change in the pod size calibration parameters the derivatives are more reliable.
\item In the package example (below), the code was added to compute the abundance series
and v-c matrix for alternative sets of matching parameters that provide a sequence of 
lower and greater abundance estimates. These will be used in the Punt-Wade assessment
paper.
\item \code{\link{gamma.pod}} was added to data to store the fitted model to prevent it from
having to be fit each time.  This was the model used to create \code{\link{gsS}}.
}
Version 1.6:
\itemize{
\item Fixed example code in ERAnalysis to use as.numeric before using cut on date difference calculations.
That code was only for data summarization and did not affect abundance estimates. The change was necessary
to accommodate a change in R. Thanks to Brian O'Gorman for notifying me of the problem.
}
}
\author{
Jeff Laake, National Marine Mammal Laboratory, Alaska Fisheries Science Center
7600 Sand Point Way NE, Seattle, WA 98115

Observation matching and linking code written by Rod Hobbs, National Marine Mammal Laboratory, Alaska Fisheries Science Center
7600 Sand Point Way NE, Seattle, WA 98115

Maintainer:  <jeff.laake@noaa.gov>

}
\references{
BUCKLAND, S. T., J. M. BREIWICK, K. L. CATTANACH and J. L. LAAKE. 1993. Estimated population size of the California gray whale. Marine Mammal Science 9: 235-249.

BUCKLAND, S. T. and J. M. BREIWICK. 2002. Estimated trends in abundance of eastern Pacific gray whales from shore counts (1967/68 to 1995/96). Journal of Cetacean Research and Management 4: 41-48.

HOBBS, R. C., D. RUGH, J. WAITE, J. M. BREIWICK and D. P. DEMASTER. 2004. Abundance of eastern North Pacific gray whales on the 1995/96 southbound migration. Journal of Cetacean Research and Management 6: 115-120.

LAAKE, J. L., D. J. RUGH, J. A. LERCZAK and S. T. BUCKLAND. 1994. Preliminary estimates of population size of gray whales from the 1992/93 and 1993/94 shore-based surveys.

LAAKE, J., A. PUNT, R. HOBBS, M. FERGUSON, D. RUGH, and J. BREIWICK. 2009. Re-evaluation of gray whale southbound surveys. NOAA Tech. Memo.

RUGH, D. J., R. C. FERRERO and M. E. DAHLHEIM. 1990. Inter-observer count discrepances in a shore-based census of gray whales (Eschrichtius robustus). Marine Mammal Science [Mar. Mamm. Sci.]. Vol. 6: 109-120.

RUGH, D. J., M. M. MUTO, R. C. HOBBS and J. A. LERCZAK. 2008. An assessment of shore-based counts of gray whales. Marine Mammal Science [Mar. Mamm. Sci.]. Vol. 24: 864-880.

RUGH, D., R. C. HOBBS, J. A. LERCZAK and J. M. BREIWICK. 2005. Estimates of abundance of the eastern North Pacific stock of gray whales (Eschrictius robustus) 1997-2002. Journal of Cetacean Research and Management 7: 1-12.

RUGH, D. J., J. M. BREIWICK, M. E. DAHLHEIM and G. C. BOUCHER. 1993. A comparison of independent, concurrent sighting records from a shore-based count of gray whales. Wildlife Society Bulletin 21: 427-437.

RUGH, D., J. BREIWICK, M. M. MUTO, R. HOBBS, K. SHELDEN, C. D'VINCENT, I. M. LAURSEN, S. REIF, S. MAHER and S. NILSON. 2008. Report of the 2006-2007 census of the eastern north Pacific stock of gray whales. NOAA Tech Report.

}
\examples{
# The recent survey data 1987 and after are stored in ERSurveyData and those data
# are processed by the ERAbund program to produce files of sightings and effort.
# The sightings files are split into Primary observer and Secondary observer sightings.
# Primary observer sightings are whales that are not travelling North and are defined by
# those when EXPERIMENT==1 (single observer) or a designated LOCATION when EXPERIMENT==2.
#  For surveys 2000/2001 and 2001/2002, the primary observer was at LOCATION=="N"
# and for all other years, LOCATION=="S".
#
# Based on the projected timing of the passage of the whale (t241) perpendicular to the 
# watch station, the sighting was either contained in the watch (on effort) or not (off effort).
# The dataframe Primary contains all of the on effort sightings and PrimaryOff contains all
# of the off-effort sightings.  The code below shows that the counts of primary
# sighting records matches the total counts of the on and off effort sightings split into
# the 2 dataframes.
#
data(PrimaryOff)
data(Primary)
data(ERSurveyData)
NorthYears=c(2000,2001)
rowSums(with(ERSurveyData[ERSurveyData$EFLAG==3&ERSurveyData$TRAVELDIR!="N"&(ERSurveyData$EXPERIMENT==1 | (ERSurveyData$EXPERIMENT==2& 
((ERSurveyData$Start.year\%in\%NorthYears & ERSurveyData$LOCATION=="N") | (!ERSurveyData$Start.year\%in\%NorthYears & ERSurveyData$LOCATION=="S")))),] ,table(Start.year,LOCATION)))
table(Primary$Start.year) + table(PrimaryOff$Start.year)
# Likewise, the secondary sightings are those with EXPERIMENT==2 but the LOCATION that
# is not designated as primary.  The following show that the counts match for ERSurveyData
# and the dataframe SecondarySightings that was created by ERAbund
data(SecondarySightings)
{rowSums(with(ERSurveyData[ERSurveyData$EFLAG==3&ERSurveyData$TRAVELDIR!="N"&
                  (ERSurveyData$EXPERIMENT==2& ( 
                  (ERSurveyData$Start.year\%in\%NorthYears & ERSurveyData$LOCATION=="S") | 
                  (!ERSurveyData$Start.year\%in\%NorthYears & (ERSurveyData$LOCATION=="N" | ERSurveyData$LOCATION=="NP"))
                  )),],table(Start.year,LOCATION)))}
table(SecondarySightings$Start.year)
{rowSums(with(ERSurveyData[ERSurveyData$EFLAG==3&ERSurveyData$TRAVELDIR!="N"&
                   ERSurveyData$VISCODE<=4&ERSurveyData$WINDFORCE<=4&
                  (ERSurveyData$EXPERIMENT==1 | 
                  (ERSurveyData$EXPERIMENT==2& ( 
                  (ERSurveyData$Start.year\%in\%NorthYears & ERSurveyData$LOCATION=="N") | 
                  (!ERSurveyData$Start.year\%in\%NorthYears & ERSurveyData$LOCATION=="S")
                  ))),],table(Start.year,LOCATION)))}
# The data in PrimarySightings are all southbound sightings for all years in which visibility and beaufort
# are less than or equal to 4. Below the counts are shown for the 2 dataframes for
# recent surveys since 1987/88.
table(Primary$Start.year[Primary$vis<=4 & Primary$beaufort<=4])
data(PrimarySightings)
table(PrimarySightings$Start.year[PrimarySightings$Start.year>=1987])
# The following code produces the series of abundance estimates used in the paper and
# some comparative naive estimates.  It will also optionally (commented out) compute
# the variance-covariance matrix for the estimates.
#
# Get sightings and effort data
data(PrimarySightings)
data(PrimaryEffort)
# The data from observer MAS in 1995 is excluded because this observer
# did not participate in the double count in that year. 
cat("\nMAS 1995/96 sightings = ",nrow(PrimarySightings[(PrimarySightings$Observer=="MAS"&PrimarySightings$Start.year==1995),]))
cat("\nMAS hours of effort in 1995/96= ",24*sum(PrimaryEffort$effort[(PrimaryEffort$Observer=="MAS"&PrimaryEffort$Start.year==1995)]))
PrimarySightings=PrimarySightings[!(PrimarySightings$Observer=="MAS"&PrimarySightings$Start.year==1995),]
PrimaryEffort=PrimaryEffort[!(PrimaryEffort$Observer=="MAS"&PrimaryEffort$Start.year==1995),]
table(PrimarySightings$Start.year)
# Define arguments used in the analysis
all.years=unique(PrimaryEffort$Start.year)
recent.years=all.years[all.years>=1987]
early.years=all.years[all.years<1987]
final.time=sapply(tapply(floor(PrimaryEffort$time),PrimaryEffort$Start.year,max),function(x) ifelse(x>90,100,90))
lower.time=rep(0,length(final.time))
fn=1.0817
se.fn=0.0338
# Effort and sightings prior to 1987 were filtered for an entire watch if vis or beaufort 
# exceeded 4 at any time during the watch.  This is done for surveys starting in 1987 with the
# Use variable which is set to FALSE for all effort records in a watch if at any time the vis or
# beaufort exceeded 4 during the watch.
# Here are the hours of effort that are excluded (FALSE) and included (TRUE) by each year
# Note that for most years <1987 there are no records with Use==FALSE because the filtered records
# were excluded at the time the dataframe was constructed. The only exception is for 1978 in which  
# one watch (5 hours) was missing a beaufort value so it was excluded.
tapply(PrimaryEffort$effort,list(PrimaryEffort$Use,PrimaryEffort$Start.year),sum)*24
# These are the number of sightings that were included/excluded based on Use
Sightings=PrimarySightings
Sightings=merge(Sightings,subset(PrimaryEffort,select=c("key","Use")))
tapply(Sightings$Start.year,list(Sightings$Use,Sightings$Start.year),length)
# Filter effort and sightings and store in dataframes Effort and Sightings
Effort=PrimaryEffort[PrimaryEffort$Use,]  
Sightings=PrimarySightings
Sightings$seq=1:nrow(Sightings)
Sightings=merge(Sightings,subset(Effort,select=c("key")))
Sightings=Sightings[order(Sightings$seq),]
# Using the filtered data, compute simple minded abundance estimates by treating the 
# sampled periods (watches) as a random sample of a period of a specified number
# of days (e.g., 4 days).  It uses raw counts with no correction for pod size or missed pods.
period=4
# compute the fraction of each period that was sampled (eg. 12 hours / (4 days *24 hrs/day))
sampled.fraction=with(Effort,
{
  Day=as.numeric(as.Date(Date)-as.Date(paste(Start.year,"-12-01",sep="")))
  tapply(effort,list(Start.year,cut(Day,seq(0,100,period))),sum)/period
})
# compute the number of whales counted in each period
whales.counted=with(Sightings,
{
  Day=as.numeric(as.Date(Date)-as.Date(paste(Start.year,"-12-01",sep="")))
  tapply(podsize,list(Start.year,cut(Day,seq(0,100,period))),sum)
})
# Compute simple minded population estimate and plot it
period.estimate=apply(whales.counted/sampled.fraction,1,sum,na.rm=TRUE)
pdf("SimpleMindedEstimates.pdf")
plot(all.years,period.estimate,xlab="Survey year",ylab="Estimated gray whale population size")
dev.off()
# Compute naive estimates of abundance for the 23 surveys.  These use the uncorrected
# counts of whales from the primary observer during watches in which neither Beaufort nor
# vis exceeded 4.  For each year a gam with a smooth over time is fitted and this is
# used to predict total abundance throughout the migration from the counts of whales
# during the sampled periods.  There is no correction for missed pods or for measurement
# error in podsize. Each fitted migration gam is plotted with the observed values and
# saved in the file NaiveMigration.pdf.
pdf("NaiveMigration.pdf")
naive.abundance.models=vector("list",23)
i=0
for (year in all.years)
{
   i=i+1
   primary=Sightings[Sightings$Start.year==year,]
   primary$Start.year=factor(primary$Start.year)
   ern=subset(Effort,subset=as.character(Start.year)==year,
         select=c("Start.year","key","begin","end","effort","time","vis","beaufort"))
   ern$Start.year=factor(ern$Start.year)
   naive.abundance.models[[i]]=estimate.abundance(spar=NULL,
     dpar=NULL,gsS=gsS,effort=ern, sightings=primary, final.time=final.time[i],lower.time=lower.time[i],
     gformula=~s(time),dformula=NULL)
}
dev.off()
Nhat.naive=sapply(naive.abundance.models,function(x) x$Total)
pdf("NaiveMigrationEstimates.pdf")
plot(all.years,Nhat.naive,xlab="Survey year",ylab="Estimated gray whale population size",pch=16)
dev.off()
# Define set of models to be evaluated for detection
models=c("podsize+Dist+Observer",
         "podsize+Dist+Observer+beaufort",
         "podsize+Dist+Observer+vis",
         "podsize+Dist+Observer+Vis")
# Create time series of estimates based on previous approach using Reilly pod size
# correction method but using 1978 data for surveys <=1987 and 1992-1994 aerial data
# for surveys >=1992
data(add.cf.reilly)
data(add.cf.laake)
Sightings$corrected.podsize[Sightings$Start.year<=1987]=reilly.cf(Sightings$podsize[Sightings$Start.year<=1987],add.cf.reilly)
Sightings$corrected.podsize[Sightings$Start.year>1987]=reilly.cf(Sightings$podsize[Sightings$Start.year>1987],add.cf.laake)
pdf("ReillyApproach.pdf")
reilly.estimates=compute.series(models,naive.abundance.models,sightings=Sightings,effort=Effort,TruePS=FALSE)
Nhat.reilly=reilly.estimates$Nhat
Nhat.reilly[1:15]=Nhat.naive[1:15]*(Nhat.reilly[16]/Nhat.naive[16])
avg.Reilly.podsize=tapply(Sightings$corrected.podsize,Sightings$Start.year,mean)
dev.off()
pdf("ReillyEstimates.pdf")
plot(all.years,Nhat.reilly,xlab="Survey year",ylab="Estimated gray whale population size",pch=16)
dev.off()
# Next compute the series of abundance estimates for the most recent 8 years by
# fitting and selecting the best detection model but not applying the pod size correction.
# From those 8 estimates and the naive estimates, compute an average ratio and 
# apply it to generate the estimates for the first 15 surveys prior to 1987.
Sightings$corrected.podsize=Sightings$podsize
abundance.estimates.nops.correction=compute.series(models,naive.abundance.models,sightings=Sightings,effort=Effort,TruePS=FALSE)
Sightings$corrected.podsize=NULL
cat("\nCorrection factor for missed pods without pod size correction\n")
print(cbind(Year=all.years,ratio=abundance.estimates.nops.correction$Nhat/Nhat.naive/fn))
# Next compute the series of abundance estimates for the most recent 8 years by
# fitting and selecting the best detection model.  From those 8 estimates and the
# naive estimates, compute an average ratio and apply it to generate the estimates
# for the first 15 surveys prior to 1987. Note with hessian=TRUE, the analysis can
# take about 30-60 minutes to complete.
pdf("Migration.pdf")
abundance.estimates=compute.series(models,naive.abundance.models,sightings=Sightings,effort=Effort,hessian=TRUE)
dev.off()
pdf("CurrentEstimates.pdf")
plot(all.years,abundance.estimates$Nhat,xlab="Survey year",ylab="Estimated gray whale population size")
dev.off()
cat("\nCorrection factor for pod size error\n")
print(cbind(Year=all.years,ratio=abundance.estimates$Nhat/abundance.estimates.nops.correction$Nhat))
#
# Print out estimates for Table 7 in Laake et al. -- pod size and detection  
  for(i in 1:8)
   print(cbind(paste(signif(abundance.estimates$det[[i]][[1]]$par,digits=3)," (",
   signif(sqrt(diag(solve(abundance.estimates$det[[i]][[1]]$model$hessian))),digits=3),")",sep="")))
#
# Next compute E(S) and it's std error for Table 7
#
 ES.var=function(par,hessian,gsS,nmax=20)
 {
    spar=par[1:2]
    vc=solve(hessian)[1:2,1:2]
    sparest=spar
    sparest[1]=spar[1]*0.999
    s1.low=expected.podsize(sparest,gsS=gsS,nmax=nmax)$ES
    sparest[1]=spar[1]*1.001
    s1.hi=expected.podsize(sparest,gsS=gsS,nmax=nmax)$ES
    sparest=spar
    sparest[2]=spar[2]*0.999
    s2.low=expected.podsize(sparest,gsS=gsS,nmax=nmax)$ES
    sparest[2]=spar[2]*1.001
    s2.hi=expected.podsize(sparest,gsS=gsS,nmax=nmax)$ES
    deriv.s=c((s1.hi-s1.low)/(.002*spar[1]),(s2.hi-s2.low)/(.002*spar[2]))
    return(list(ES=expected.podsize(spar,gsS=gsS,nmax=nmax)$ES,se=sqrt(t(deriv.s)\%*\%vc\%*\%deriv.s)))
}
  data(gsS)
  for(i in 1:8)
    print( ES.var(abundance.estimates$detection.models[[i]][[1]]$model$par,abundance.estimates$detection.models[[i]][[1]]$model$hessian,gsS))
# NOT RUN: Compute the var-cov matrix of the estimates.  These are not run by default
# because it takes a very long time to run. To run copy code from documentation into R.
\dontrun{
#
# Load the fitted mixed-effects gamma model with a random pod effect -- this was the most parsimonious model
 data(gamma.pod)
#
 ps.results=gamma.pod
 ps.results$vc=solve(ps.results$hessian)
 abundance.vc=compute.series.var(abundance.estimates,naive.abundance.models,ps.results=ps.results,delta=c(0.001,0.01))
 write.table(log(abundance.estimates$Nhat),"nhat.txt")
 lnvc=log(1+abundance.vc$vc/outer(abundance.estimates$Nhat,abundance.estimates$Nhat,"*"))
 write.table(lnvc,"vc.txt")
# The following produces Fig 4 in Laake et al.
 pdf("EstimatesWithCL.pdf")
 gw.abund <-
   structure(list(year = c(1967, 1968, 1969, 1970, 1971, 1972, 1973,
   1974, 1975, 1976, 1977, 1978, 1979, 1984, 1985, 1987, 1992, 1993,
   1995, 1997, 2000, 2001, 2006), N.old = c(13776, 12869, 13431, 11416,
   10406, 16098, 15960, 13812, 15481, 16317, 17996, 13971, 17447,
   22862, 21444, 22250, 18844, 24638, 24065, 29758, 19448, 18178,
   20110), se.N.old = c(1082, 708, 758, 590, 614, 834, 872, 781,
   930, 818, 1249, 753, 984, 1379, 1120, 1115, 1190, 1475, 1393,
   3122, 1882, 1780, 1766)), .Names = c("year","N.old",
   "se.N.old"), row.names = c(NA, -23L), class = "data.frame")
 conf.int=function(abundance, CV, alpha=0.05, digits=2, prt=FALSE)
 {
# Computes confidence intervals based on lognormal distr.
# JMB / NMML / 11 Sep 2008

   if (alpha <0 || alpha > .999) stop("alpha must be in (0,1)")
   z = round(abs(qnorm(alpha/2)),2)
   if (prt) cat("N:",abundance,"  cv:",CV,"  alpha:",alpha,"  z:",z,"\n")
   C <- exp(z * sqrt(log(1 + CV^2)))
   SL <- round(abundance/C,digits)
   SU <- round(abundance * C,digits)
   data.frame(SL,SU)
  }
  gwplot=function (N,se,year,bar.wid=.5,p=NULL)
  {
#   Function written by JMB with mods by JLL
    require("ggplot2")
# Plot gray whale abund. data using ggplot2 package
    if (length(N) != length(se) | length(se) != length(year)) stop("error in
    vector lengths")
    cv <- se/N 
    yy <- conf.int(N,cv)
    y.max <- yy$SU ; y.min <- yy$SL
    yy <- cbind(y.min,y.max)
    gwd <- data.frame(year,N) 
    gwd <- cbind(gwd,yy)
    names(gwd) <- c("Year","Abundance","y.min","y.max")
    if(is.null(p))
    {
      limits <- aes(ymax = y.max, ymin = y.min)   
      p <- ggplot(gwd,aes(x=Year,y=Abundance))
      p <- p + geom_point(colour="darkblue",data=gwd) +
      geom_errorbar(limits,colour="blue",width=bar.wid,data=gwd)
    }
    else
    {
      limits <- aes(ymax = y.max, ymin = y.min)   
      p <- p + geom_point(colour="black",data=gwd) +
      geom_errorbar(limits,colour="black",linetype=2,width=bar.wid,data=gwd)
    }
    p
  }
  xx=gwplot(abundance.estimates$Nhat,abundance.vc$se,gw.abund$year)
  gwplot(gw.abund$N.old,gw.abund$se.N.old,gw.abund$year-.5,p=xx)
 dev.off()
 # Alternate runs using different matching weights
  abundance.estimates.lo=compute.series(models,naive.abundance.models,sightings=Sightings,effort=Effort,hessian=TRUE,twt=0.11,dwt=3.02)
  write.table(log(abundance.estimates.lo$Nhat),"nhat_lo.txt")
  abundance.estimates.hi=compute.series(models,naive.abundance.models,sightings=Sightings,effort=Effort,hessian=TRUE,twt=0.27,dwt=5.06)
  write.table(log(abundance.estimates.hi$Nhat),"nhat_hi.txt")
  abundance.vc.lo=compute.series.var(abundance.estimates.lo,naive.abundance.models,ps.results=ps.results,twt=0.11,dwt=3.02,delta=c(0.001,0.01))
  lnvc=log(1+abundance.vc.lo$vc/outer(abundance.estimates.lo$Nhat,abundance.estimates.lo$Nhat,"*"))
  write.table(lnvc,"vc_lo.txt")
  abundance.vc.hi=compute.series.var(abundance.estimates.hi,naive.abundance.models,ps.results=ps.results,twt=0.27,dwt=5.06,delta=c(0.001,0.02))
  lnvc=log(1+abundance.vc.hi$vc/outer(abundance.estimates.hi$Nhat,abundance.estimates.hi$Nhat,"*"))
  write.table(lnvc,"vc_hi.txt")
}
# Plot for each year average observed pod size and previously used corrected pod size 
  data(gsS)
  pdf("PodsizeComparisons.pdf")
  avg.podsize=tapply(Sightings$podsize,Sightings$Start.year,mean)
  plot(all.years, avg.podsize, ylim=c(1,3), xlab="Survey Year (1 Dec)",ylab="Average gray whale pod size",pch=1)
  points(all.years, avg.Reilly.podsize, pch=2)
  legend(1970,1.25,pch=c(1,2),legend=c("Observed average","Reilly correction"))
# Next plot the averages for the last 8 years after linking pods, Reilly corrected value and then the expected value
# from the fitted gamma distribution
  plot(recent.years, abundance.estimates$summary.df$MeanPS, ylim=c(1,3), xlab="Survey Year (1 Dec)",ylab="Average gray whale pod size",pch=1)
  expected.ps=sapply(abundance.estimates$detection.models,function(x) expected.podsize(x[[1]]$par[1:2],gsS,20)$ES)
  points(recent.years,expected.ps, pch=16)
  points(recent.years,avg.Reilly.podsize[16:23], pch=2)
  legend(1990,1.35,pch=c(1,2,16),legend=c("Observed average","Reilly correction","Expected pod size"))
  dev.off()
# Plot distributions of pod size in calibration data and for true pod size from that year
# This is Figure 5 in Laake et al.
data(PodsizeCalibrationTable)
obs.1992=with(PodsizeCalibrationTable[substr(PodsizeCalibrationTable$key,1,11)=="Aerial_1992",],table(cut(True,c(1,2,3,4,20),right=FALSE)))
obs.1993=with(PodsizeCalibrationTable[substr(PodsizeCalibrationTable$key,1,11)=="Aerial_1993",],table(cut(True,c(1,2,3,4,20),right=FALSE)))
obs.1997=with(PodsizeCalibrationTable[substr(PodsizeCalibrationTable$key,1,8)=="Tracking",],table(cut(True,c(1,2,3,4,20),right=FALSE)))
obs.1992=obs.1992/sum(obs.1992)
obs.1993=obs.1993/sum(obs.1993)
obs.1997=obs.1997/sum(obs.1997)
ps.1992=gammad(c(-.073,-.3474),20)
ps.1993=gammad(c(-.07,-.474),20)
ps.1997=gammad(c(-.598,-.674),20)
ps.1992=c(ps.1992[1:3],sum(ps.1992[4:20]))
ps.1993=c(ps.1993[1:3],sum(ps.1993[4:20]))
ps.1997=c(ps.1997[1:3],sum(ps.1997[4:20]))
par(mfrow=c(3,1))
barplot(rbind(obs.1992,ps.1992),beside=TRUE,space=c(0,.1),main="1992-93",ylab="Proportion")
barplot(rbind(obs.1993,ps.1993),beside=TRUE,space=c(0,.1),main="1993-94",ylab="Proportion")
barplot(rbind(obs.1997,ps.1997),beside=TRUE,space=c(0,.1),main="1997-98",ylab="Proportion",xlab="True pod size",legend.text=c("Calibration","Estimated"),args.legend=list(x=8,y=.8))
# Code to construct Figure 2
par(mfrow=c(2,2))
data(PodsizeCalibrationTable)
psdf=PodsizeCalibrationTable
ps1=as.matrix(psdf[psdf$True==1,3:22])
ps2=as.matrix(psdf[psdf$True==2,3:22])
ps3=as.matrix(psdf[psdf$True==3,3:22])
psplus=as.matrix(psdf[psdf$True>=4,3:22])
gpar=gamma.pod$par
ps.results=gamma.pod
expect.ps=create.gsS(ps.results,True=psdf$True[psdf$True>=4][1:21])
expect.ps=colMeans(expect.ps)
barplot(rbind(gsS[1,],colMeans(ps1/rowSums(ps1))),beside=TRUE,main="True size = 1",legend.text=c("exp","obs"))
barplot(rbind(gsS[2,],colMeans(ps2/rowSums(ps2))),beside=TRUE,main="True size = 2",legend.text=c("exp","obs"))
barplot(rbind(gsS[3,],colMeans(ps3/rowSums(ps3))),beside=TRUE,main="True size = 3",legend.text=c("exp","obs"))
barplot(rbind(expect.ps,colMeans(psplus/rowSums(psplus))),beside=TRUE,main=paste("True size = 4+",sep=""),legend.text=c("exp","obs"))
}


