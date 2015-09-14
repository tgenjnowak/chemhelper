\name{merge_scan_events}
\alias{merge_scan_events}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
merge even numbered scans with preceding odd numbered scan
}
\description{
merge even numbered scans with preceding odd numbered scan. This is used if for example the m/z range was split between different traces ("functions" in the lingo of some vendors).
}
\usage{
merge_scan_events(infile,outfile,oddrange,evenrange)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{infile}{Path to any file readable by \code{\link{xcmsRaw}}.}
  \item{outfile}{Filename of output mzdata file.}
  \item{oddrange}{numeric vector of length 2 giving the min and max m/z values to be extracted from the odd numbered scans.}
  \item{evenrange}{numeric vector of length 2 giving the min and max m/z values to be extracted from the even numbered scans.}

}
\value{None}
\author{
Jan Stanstrup, \email{stanstrup@gmail.com}
}
