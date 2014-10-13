\name{smile2inchi}
\alias{smile2inchi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert SMILES to InChIs
}
\description{
Convert SMILES to InChIs. 
Calls obabel locally. Hence obabel need to be available in path.
}
\usage{
smile2inchi(smile,verbose=F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{smile}{Character vector SMILES.}
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
