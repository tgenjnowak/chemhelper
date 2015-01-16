\name{orbifilter,xcmsRaw.orbifilter}
\alias{orbifilter}
\alias{xcmsRaw.orbifilter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Remove satellite/shoulder peaks (orbitrap artifacts)
}
\description{
The orbifilter function takes a matrix with columns "mz" and "intensity" (single spectrum) and xcmsRaw.orbifilter takes an xcmsRaw object and tries to remove artifact peaks caused by fourier transformation and centroiding on orbitrap instruments.

In each scan the functions looks for the largest peak. It then looks in the range windows_width/2 around this peak for peaks with a relative intensity below max_rel_int.
 
Optionally it can try to detect if a peak is really an isotope (could have multiple charges and thus have mass close to the real peak). This is enabled with keep_isotopes. It then looks for the equivalent difference that would be associated with isotopes for molecules with a charge between 1 and max_charge. isotope_mz_tol adjusts the requirement for the accuracy of the mass difference to the isotope.

The process is continues untill all peaks have either been marked for removal or all peaks have been assessed for neighboring artifacts.

}
\usage{
orbifilter(x,windows_width=0.1,max_rel_int = 0.2,keep_isotopes=TRUE,max_charge=5,isotope_mz_tol = 0.005)
xcmsRaw.orbifilter(xraw,windows_width=0.1,max_rel_int = 0.2,keep_isotopes=TRUE,max_charge=5,isotope_mz_tol = 0.005)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Matrix or data.frame with columns "mz" and "intensity".}
  \item{xraw}{\code{\link{xcmsRaw-class}} object.}
  
  \item{windows_width}{Window for locating satellite/shoulder peaks.}
  \item{max_rel_int}{Max relative intensity for a peak to be identified as a atellite/shoulder peak.}
  \item{keep_isotopes}{ Keep masses if they could be isotopes. Also multicharged. }
  \item{max_charge}{Max charge for locating possible isotope peaks.}
  \item{isotope_mz_tol}{Adjusts the requirement for the accuracy of the mass difference to the possible isotope.}
  
  
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{matrix/data.frame or \code{\link{xcmsRaw-class}}.}
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
