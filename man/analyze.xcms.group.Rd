\name{analyze.xcms.group}
\alias{analyze.xcms.group}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot and analyze xcms peak groups in rt / m/z space
}
\description{
Plots all peaks (in all samples) in a given rt / m/z area. The peak group boundaries are draw as rectangles. A large red spot markes the point of interest as defined in the input to the function. This can be used to debug peak groups that are erroneously split or merged. 
}
\usage{
analyze.xcms.group(xset_grouped,mz,rt,
                   rt_tol_sample=60,
                   mz_tol_sample=0.01,
                   rt_tol_group=30,
                   mz_tol_group=0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xset_grouped}{\code{\link{xcmsSet-class}} to analyze.}
  \item{mz}{Approximate m/z of the peak(s) of interest.}
  \item{rt}{Approximate retention time (seconds) of the peak(s) of interest.}
  \item{rt_tol_sample}{Retention time (seconds) tolerances for finding related peaks across all samples.}
  \item{mz_tol_sample}{m/z tolerances for finding related peaks across all samples.}
  \item{rt_tol_group}{Retention time (seconds) tolerances for finding related peaks groups.}
  \item{mz_tol_group}{m/z tolerances for finding related peaks groups.}
}

\value{A matrix giving the peak groups in the selected rt / m/z area.}

\author{
Jan Stanstrup, \email{stanstrup@gmail.com}
}
\examples{
library(xcms)
library(faahKO)
filepath <- system.file("cdf", package = "faahKO")
xsg <- group(faahko)
xsg <- fillPeaks(xsg)
analyze.xcms.group(xsg,mz=256.1500,rt=3451.305,
                   rt_tol_sample=300,
                   mz_tol_sample=0.01,
                   rt_tol_group=300,
                   mz_tol_group=0.05)
}
