\name{Maintenance}
\Rdversion{1.1}
\alias{Maintenance}
\title{Creation and Maintenance of the ERAnalysis Data}
\description{Details on the creation and maintenance of the dataframes in the ERAnalysis package and the links with the ERAbund program
and the ACCESS database.}

\details{
The southbound gray whale surveys have been conducted from 1967-2006 in California.  The surveys were conducted at Yankee Point
from 1967-1973 and at Granite Canyon since 1974.  The formats for the survey data changed over time but a relatively
common format was used for the years 1967-1985 and a much improved format since 1987.  Throughout the documentation the term earlier
surveys refers to the period 1967-1985 and the later period is called the more recent surveys.  The data for the earlier surveys
were originally maintained on punch cards in a variety of formats.  These data in their original various formats are maintained
in the NMML CAEP directory Gray_Wh/1 Master database gray whales/Historic data in a zip file (ERDataReFormat1967-1985.zip).  The
content of the files have been changed in places where keypunch errors or other anomalies were discovered.  All changes were documented
in the files on the right hand side of each line that was modified. These data were structured based on different card types.  A 
card type 2 was an effort record and a card type 3 was a sighting record and a card type 4 was a comment. Also, contained in the 
zip file is an R script (ReFormatOldData.R) that was used to read in the various formats and create a common effort file (card 2) and a common sighting file (card 3)
for these early surveys.  See the documentation in the R script for details.  These 2 files were transferred to the gray whale
Access database in Gray_Wh/1 Master database gray whales/GrayWhaleSurveyData.accdb and are the tables EffortPre1987 and ERSightingsPre1987.
Both of these tables are contained in dataframes in this package named \code{\link{EarlyEffort}} and \code{\link{EarlySightings}} respectively.  

The format for the more recent surveys was radically different which made it difficult and not really desirable to merge the
original data from the 2 periods.  The format for the more recent surveys is also based on card-types of a sort called an 
event flag (EFLAG) but unlike the earlier survey data all of the records in the data have a common format rather than having a different format for each card-type (EFLAG).  
Various methods of recording and storing data were used from 1987 onwards but the format remained relatively unchanged such that it
was relatively easy to merge them into a common table in the Access database named \code{AllRecentData} and which is called \code{\link{ERSurveyData}} in
this package.  

The function \code{\link{store.ERdata}} extracts the above data and saves them in the data directory of the ERAnalysis package.  If
the data in the ACCESS database are changed, then \code{\link{store.ERdata}} should be called and the package should be rebuilt to include
the changes in the data that accompany the package.  In
addition to the above dataframes, the function also saves the \code{\link{Observer}} table and the \code{\link{PodsizeCalibrationData}} and 
\code{\link{PodsizeCalibrationTable}} data tables.  These latter data are not in the ACCESS database but are in the PSCalibAll.txt file in
the package data directory.  If for any reason the platform (observer) eye heights were changed, they should be changed in the Height.txt file 
which is also in the package data directory. 

Starting with the 2000/2001 survey a program written by Jim Lerczak called ERAbund was used to record and error check data
and to process the data into various output files for effort and sightings.  All of the recent survey data from 1987 onward can and has been formatted for the ERAbund program 
to provide a uniform creation of the files.  ERAbund uses a DBase format(.dbf) for the input files and it creates text output files.  
It uses 2 input files with the extension .in that define file names and various parameters that are used for data processing 
and creation of the output files. Because these input files have specific names and cannot be specified to the program, it is necessary to store the various surveys into
separate sub-directories with the input .dbf file and the appropriate .in files for that survey year.  To automate this to some degree, the 
package assumes the sub-directories are in a common directory and are named with the format ERSWYYYY where YYYY is the last 2 digits of the
years that the survey spans (eg., 8788 is the 1987/1988 survey).  Currently, there are 8 such sub-directories for the 1987/88,1992/93,1993/94,1995/96,
1997/98, 2000/01, 2001/02, and 2005/06 surveys.  

The function \code{\link{CreateERAbundFiles}} will create all of the .dbf files from the \code{\link{ERSurveyData}} and will put them into 
the appropriately named subdirectories of the directory that you can specify and which defaults to 
c:/gw.  Each sub-directory contains the .in files needed for that survey year and these should not be changed without a good reason.  Each
assumes that all beauforts and visibility conditions are valid because this filtering can be done with the code contained in this package.  
They only differ in slight ways like platform height which has varied depending on the year and the primary observer designation that was 
the observer in the south location for most years except for the 2000/01 and 2001/02 surveys. To get the analysis output files, 
it is necessary to run the ERAbund program in each sub-directory which is most easily done fromn a command (DOS) window after adding 
the directory for the ERAbund.exe to your path.  ERAbund is a Visual-Basic program that presents a screen and menu.  You only need to 
choose the analysis menu option and then close all of the windows once it has completed.  
It creates a series of output files that are documented in the ERAbund User's Manual.  The files used with this package are
named as: 1) Primary observer sighting file: erprimyyyy.dat, 2) Secondary observer sighting file: ersecyyyy.dat (not surrently used by the code in this package), 3) Match file for double-observer data: 
ermatchspyyyy.dat,and 4) Effort data for fitting the migration curve: ernormyyyy.dat.  The function \code{\link{ReadERAbundFiles}} reads in
all of these files from the sub-directories for the specified years and creates a dataframe for each file. To get the secondary effort files, 
a separate Secondary sub-directory was created and under it are the individual ERSWYYY subdirectories but the .in files are set such that the 
secondary location is treated as the primary location to get all of the secondary effort.

If for the data in the ACCESS database are modified for surveys in 1987 or afterwards, the the following steps should be
taken:

1) Run \code{\link{CreateERAbundFiles}}; note that a warning message will be issued and can be ignored.  As long as the directory structure is
in c:/gw and the ACCESS database is in the current R working directory then no arguments are needed.

2) Run the analysis menu option of ERAbund in each of the primary and secondary subdirectories to create the output files.

3) Run \code{\link{store.ERAbundData}} to extract the contents of the output files to create the following package dataframes:
\code{\link{Primary}}, \code{\link{PrimaryOff}}, \code{\link{PrimarySightings}}, \code{\link{SecondarySightings}}, \code{\link{Match}},
\code{\link{PrimaryEffort}}, \code{\link{AllPrimaryEffort}}, \code{\link{SecondaryEffort}}, \code{\link{AllSecondaryEffort}}.

4) Rebuild the package to incorporate the revised dataframes into the package.
 
To make the analysis process easier, the effort and sightings data from the early and recent surveys are merged into common dataframes using a minimal
amount of common data fields that are available for all surveys years.  These dataframes are \code{\link{PrimaryEffort}} and
\code{\link{PrimarySightings}} which are the selected effort and sighting records for the primary observer which do not exceed specified maximum
values for the visibility and Beaufort codes.  This data selection and formatting is done by the function \code{\link{process.ERdata}} and it is
called when \code{\link{store.ERAbundData}} is executed.  In addition, it also creates \code{\link{AllPrimaryEffort}} which is the merged
effort without beaufort/vis conditions.  This latter file and \code{\link{AllSecondaryEffort}} are only used in the matching process
by \code{\link{create.match}}.  

Note that while the dataframe \code{\link{Match}} created by ERAbund is still maintained in the
package, those data are not used in the analysis because there was no way to automate changes in matching parameters and a few minor errors
were discovered in the ERAbund matching code (e.g., unless both observers started at the exact same time the data were not used). 
Instead, \code{\link{create.match}} and accompanying functions were created to build the match data.  It provides much more flexibility in
the type of criterion and in the weights used in the criterion.  

The only other potential data changes would be with the pod size calibration data.  The text file PsCalibAll.txt contains
the raw data and it is reformed in \code{\link{PodsizeCalibrationData}} and \code{\link{PodsizeCalibrationTable}}.  Also, the function
\code{\link{create.podsize.calibration.matrix}} fits models to those data and creates the dataframe \code{\link{gsS}} which is contained
in the package.  In addition, that function also creates \code{\link{add.cf.all}},\code{\link{add.cf.reilly}},\code{\link{add.cf.laake}} which 
are the vectors of additive podsize correction factors.  If any of the pod size calibration data were changed, \code{\link{create.podsize.calibration.matrix}}
would also have to be re-run in addition to \code{\link{store.ERdata}} and the package would have to be re-built to incorporate those
changes.
}
\author{Jeff Laake}



