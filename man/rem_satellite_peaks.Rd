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
}
\usage{
rem_satellite_peaks(files,outdir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{files}{Paths to files to process.}
  \item{outdir}{Output folder path. Files are written to mzML format by openMS.}
  
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{None}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Jan Stanstrup, \email{stanstrup@gmail.com}
}
\note{
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
