\name{xsAnnotate_stats}
\alias{xsAnnotate_stats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
xsAnnotate statistics.
}
\description{
Give some statistics from an xsAnnotate object annotated by CAMERA.
}
\usage{
xsAnnotate_stats(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{\link{xsAnnotate}} to analyze.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A matrix giving the following statistics: \cr
'# Feature' \cr
'# Feature groups' \cr
'# Features groups with several features' \cr
'# Features in feature groups with several features' \cr
'# Features annotated as isotopes' \cr
'# Feature with adducts/fragments annotation' \cr
'# Feature groups with fragment/adduct annotation' \cr
'# Annotation groups' \cr
'Mean number of annotation groups per feature group with annotations' \cr
'Median number of annotation groups per feature group with annotations'
}
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
 library(CAMERA)
 file <- system.file('mzdata/MM14.mzdata', package = "CAMERA")
 xs   <- xcmsSet(file, method="centWave", ppm=30, peakwidth=c(5,10))
 an   <- xsAnnotate(xs)
 an   <- groupFWHM(an)
 an   <- findIsotopes(an)  # optional but recommended.
 an   <- groupCorr(an) # optional but very recommended step
 an   <- findAdducts(an,polarity="positive")
 xsAnnotate_stats(an)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
