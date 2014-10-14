\name{inchi.keep.cont}
\alias{inchi.keep.cont}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Keep only largest continues part of a molecule InChI (AKA remove salts)
}
\description{
Keep only largest continues part of a molecule InChI. This is most useful to neutralize salts.
Calls obabel locally. Hence obabel need to be available in path.
}
\usage{
inchi.keep.cont(inchi,verbose=F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inchi}{Character vector InChIs.}
  \item{verbose}{Print additional information from obabel.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{A character vector of InChIs.}
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
