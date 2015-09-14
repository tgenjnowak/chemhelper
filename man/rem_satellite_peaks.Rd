\name{rem_satellite_peaks}
\alias{rem_satellite_peaks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Remove satellite/shoulder peaks (orbitrap artifacts)
}
\description{
Remove satellite/shoulder peaks (orbitrap artifacts).
openMS need to be installed and SpectraFilterWindowMower available in path.
The SpectraFilterWindowMower program is used with the parameters: -algorithm:windowsize 0.3 and -algorithm:peakcount 1.

This function is now superseeded by \code{\link{xcmsRaw.orbifilter}} that works better and does not require openMS.
}
\usage{
rem_satellite_peaks(files,outdir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{files}{Paths to files to process.}
  \item{outdir}{Output folder path. Files are written to mzML format by openMS.}
}

\value{None}

\author{
Jan Stanstrup, \email{stanstrup@gmail.com}
}
