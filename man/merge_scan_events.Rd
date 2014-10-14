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
