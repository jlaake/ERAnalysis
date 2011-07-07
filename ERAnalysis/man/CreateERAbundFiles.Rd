\name{CreateERAbundFiles}
\Rdversion{1.1}
\alias{CreateERAbundFiles}
\title{Creates dbf files for ERAbund data entry/analysis program}
\description{
Creates a dbf file for each specified year by extracting results from the
gray whale Acess database.  The function assumes that the correct directory structure
already exists with a separate subdirectory for each survey year. Each sub-directory contains
the appropriate .in files for that survey year which set parameters for the ERAbund program.}
\usage{
CreateERAbundFiles(directory="c:/gw",DBDirectory="",years=c(87,92,93,95,97,0,1,6))
}
\arguments{
  \item{directory}{directory containing sub-directories for ERAbund dbf analysis with
                   one sub-directory for each survey}
  \item{DBDirectory}{directory containing Access database; default assumes it is in
                      the workspace directory}
  \item{years}{vector of last 2 digits of survey year; 87 is 1987/88 survey and 0 is 2000/2001 survey;
               the subdirectories are named ERSWyyyy where yyyy is the last 2 digits of
               the consecutive years like 8788}
}
\details{The dbf files are created as ERSWyyyy.dbf and stored in the sub-directory with
the same name.  The .in file should contain the same naming convention so ERAbund will
refer to the correct file.

Note that you will get the following warning message with this function is used but it can be ignored:

Warning message:
In max(nchar(x[!is.na(x)], "b")) :
  no non-missing arguments to max; returning -Inf

}
\value{None}
\author{Jeff Laake}
\seealso{\code{\link{ReadERAbundFiles}}}
